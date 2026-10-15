#!/usr/bin/env python3
"""Step 2: single point comparison on the unrelaxed canonical clean slab,
open boundaries against periodic boundaries.

EXECUTION STATUS: NOT RUN. This script has been written and reviewed but has
not been executed. Every number it prints will be produced on your machine.

Reports, for each boundary condition that the chosen method supports:
  - total energy
  - maximum force on free (unfrozen) atoms, using the frozen optimiser
    constraint FixAtoms(z < z_freeze)
  - force magnitudes on the named edge O atoms against interior O atoms

Atom labels follow the convention confirmed from inputs/state_set.json: the
site names o21, o45, o69, o93 are zero based ASE atom indices, not O ordinals.
The edge list below uses the same convention.

Usage:
    python scripts/singlepoint_edge_vs_interior.py \
        --slab inputs/slab_clean_2x2_canonical.in \
        --method GFN2-xTB \
        --out outputs/step2_edge_vs_interior.json

Two different failure modes are expected here, and they need different
handling:

  xtb-python under PBC raises a catchable Python exception, ultimately
  "scf: Multipoles not available with PBC".

  tblite with GFN2-xTB segfaults on at least macOS arm64, even for a
  non-periodic water molecule. A SIGSEGV is an OS-level crash, not an
  exception, so no try/except can trap it and an in-process loop would lose
  every row measured before it, including the xtb-python rows that
  succeeded.

Each trial therefore runs in its own subprocess, exactly as
pbc_feasibility_probe.py does, and the output JSON is rewritten after every
row. Both failure modes end up as recorded results, and a crash costs only
the row it happened on.
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

import numpy as np
from ase.constraints import FixAtoms
from ase.io import read

# The eight edge O atoms named in the audit, as zero based ASE indices.
EDGE_O = [65, 70, 82, 89, 14, 16, 56, 63]

PBC_CHOICES = {
    "FFF": (False, False, False),
    "TTF": (True, True, False),
    "TTT": (True, True, True),
}


def exception_chain(exc: BaseException) -> list:
    """Every message in a chained exception, outermost first. ASE hides the
    calculator library's own diagnosis one level down in __context__, so
    str(exc) alone returns a generic wrapper message."""
    msgs, seen, cur = [], set(), exc
    while cur is not None and id(cur) not in seen:
        seen.add(id(cur))
        msgs.append(f"{type(cur).__name__}: {cur}".strip())
        cur = cur.__cause__ or cur.__context__
    return msgs


def make_calculator(kind: str, method: str):
    if kind == "xtb-python":
        from xtb.ase.calculator import XTB
        return XTB(method=method)
    if kind == "tblite":
        from tblite.ase import TBLite
        return TBLite(method=method)
    raise ValueError(f"unknown calculator kind: {kind}")


def edge_geometry(atoms) -> dict:
    """Distance of every atom to the nearest open boundary plane in a and b.

    With pbc off, the a and b cell faces are free surfaces. For an orthorhombic
    cell the relevant quantity is simply how close an atom sits to x = 0 or
    x = a, and to y = 0 or y = b. An atom within roughly one bond length of a
    face has lost neighbours that periodicity would have supplied.
    """
    cell = atoms.get_cell()
    a_len, b_len = float(cell.lengths()[0]), float(cell.lengths()[1])
    x, y = atoms.positions[:, 0], atoms.positions[:, 1]
    d_a = np.minimum(np.abs(x - 0.0), np.abs(a_len - x))
    d_b = np.minimum(np.abs(y - 0.0), np.abs(b_len - y))
    d_edge = np.minimum(d_a, d_b)
    return {
        "a_A": a_len,
        "b_A": b_len,
        # The deepest possible interior point of an orthorhombic cell is its
        # centre, so this is the largest edge distance the cell can offer.
        "max_possible_edge_distance_A": min(a_len, b_len) / 2.0,
        "d_edge_A": d_edge,
        "d_a_A": d_a,
        "d_b_A": d_b,
    }


def single_point(slab: str, kind: str, method: str, pbc_key: str,
                 z_freeze: float) -> dict:
    atoms = read(slab, format="espresso-in")
    atoms.set_pbc(PBC_CHOICES[pbc_key])
    mask = atoms.positions[:, 2] < z_freeze
    atoms.set_constraint(FixAtoms(mask=mask))

    rec = {"calculator": kind, "method": method, "pbc": pbc_key,
           "n_frozen": int(mask.sum()), "n_free": int((~mask).sum())}
    try:
        atoms.calc = make_calculator(kind, method)
        rec["energy_eV"] = float(atoms.get_potential_energy())
        # Raw forces, before the constraint zeroes the frozen rows, so that the
        # force actually felt by a frozen edge atom is still visible.
        f = np.asarray(atoms.get_forces(apply_constraint=False), dtype=float)
        fn = np.linalg.norm(f, axis=1)

        geo = edge_geometry(atoms)
        syms = np.array(atoms.get_chemical_symbols())
        is_O = syms == "O"
        near_edge = geo["d_edge_A"] < 1.5
        interior_O = is_O & ~near_edge
        edge_O_named = np.zeros(len(atoms), dtype=bool)
        edge_O_named[EDGE_O] = True

        rec.update({
            "status": "ok",
            "fmax_free_atoms_eV_per_A": float(fn[~mask].max()),
            "fmax_all_atoms_eV_per_A": float(fn.max()),
            "named_edge_O": {
                str(i): {
                    "symbol": syms[i],
                    "d_edge_A": float(geo["d_edge_A"][i]),
                    "frozen": bool(mask[i]),
                    "force_eV_per_A": float(fn[i]),
                } for i in EDGE_O
            },
            "named_edge_O_mean_force_eV_per_A": float(fn[edge_O_named].mean()),
            "interior_O_count": int(interior_O.sum()),
            "interior_O_mean_force_eV_per_A": (
                float(fn[interior_O].mean()) if interior_O.any() else None),
            "interior_O_max_force_eV_per_A": (
                float(fn[interior_O].max()) if interior_O.any() else None),
            "geometry": {
                "a_A": geo["a_A"], "b_A": geo["b_A"],
                "max_possible_edge_distance_A": geo["max_possible_edge_distance_A"],
                "n_atoms_within_1.5A_of_edge": int(near_edge.sum()),
                "n_atoms_within_2.3A_of_edge": int((geo["d_edge_A"] < 2.3).sum()),
            },
        })
    except BaseException as exc:  # noqa: BLE001
        chain = exception_chain(exc)
        rec.update({
            "status": "failed",
            "error_type": type(exc).__name__,
            "error": str(exc),
            "error_chain": chain,
            "error_root": chain[-1] if chain else str(exc),
            "traceback_tail": traceback.format_exc().strip().splitlines()[-5:],
        })
    return rec


def run_trial_subprocess(slab: str, config: str, kind: str, method: str,
                         pbc_key: str) -> dict:
    """Run one trial in a separate interpreter.

    tblite's GFN2 path aborts the process rather than raising, so this is the
    only way the surviving rows reach the output file. returncode -11 is
    SIGSEGV, -6 is SIGABRT.
    """
    import subprocess
    proc = subprocess.run(
        [sys.executable, __file__, "--single", "--slab", slab,
         "--config", config, "--kind", kind, "--pbc-key", pbc_key,
         "--method", method],
        capture_output=True, text=True, timeout=3600,
    )
    for line in proc.stdout.splitlines():
        if line.startswith("@@RESULT@@"):
            return json.loads(line[len("@@RESULT@@"):])
    rc = proc.returncode
    signal_name = ""
    if rc < 0:
        import signal as _sig
        try:
            signal_name = _sig.Signals(-rc).name
        except ValueError:
            signal_name = f"signal {-rc}"
    return {
        "calculator": kind, "method": method, "pbc": pbc_key,
        "status": "crashed", "returncode": rc, "signal": signal_name,
        "stderr_tail": proc.stderr.strip().splitlines()[-6:],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--slab", default="inputs/slab_clean_2x2_canonical.in")
    ap.add_argument("--config", default="inputs/frozen_optimiser_config.json")
    ap.add_argument("--method", default="GFN2-xTB")
    ap.add_argument("--calculators", nargs="+", default=["xtb-python", "tblite"],
                    choices=["xtb-python", "tblite"])
    ap.add_argument("--out", default="outputs/step2_edge_vs_interior.json")
    ap.add_argument("--single", action="store_true",
                    help="Internal: run one trial and print a @@RESULT@@ line.")
    ap.add_argument("--kind", default="xtb-python")
    ap.add_argument("--pbc-key", default="FFF", choices=sorted(PBC_CHOICES))
    args = ap.parse_args()

    cfg = json.loads(Path(args.config).read_text())
    z_freeze = float(cfg["constraint"]["z_freeze_A"])

    if args.single:
        rec = single_point(args.slab, args.kind, args.method, args.pbc_key,
                           z_freeze)
        print("@@RESULT@@" + json.dumps(rec), flush=True)
        return

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    report = {"step": 2, "slab": args.slab, "method": args.method,
              "z_freeze_A": z_freeze, "rows": []}
    for kind in args.calculators:
        for pbc_key in ("FFF", "TTF", "TTT"):
            rec = run_trial_subprocess(args.slab, args.config, kind,
                                       args.method, pbc_key)
            report["rows"].append(rec)
            # Rewritten every row, so a crash never costs the earlier results.
            out.write_text(json.dumps(report, indent=2))
            if rec["status"] == "ok":
                print(f"[ok]      {kind:<11} {pbc_key}  "
                      f"E={rec['energy_eV']:.6f} eV  "
                      f"fmax_free={rec['fmax_free_atoms_eV_per_A']:.4f}  "
                      f"edgeO={rec['named_edge_O_mean_force_eV_per_A']:.4f}  "
                      f"intO={rec['interior_O_mean_force_eV_per_A']}", flush=True)
            elif rec["status"] == "crashed":
                print(f"[crashed] {kind:<11} {pbc_key}  {rec['signal']}",
                      flush=True)
            else:
                print(f"[failed]  {kind:<11} {pbc_key}  {rec['error_root']}",
                      flush=True)

    print(f"[ok] wrote {out}")


if __name__ == "__main__":
    main()
