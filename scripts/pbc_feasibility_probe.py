#!/usr/bin/env python3
"""Step 1 feasibility probe: can GFN2-xTB be evaluated on the canonical
IrO2(110) slab with periodic boundary conditions, via xtb-python or tblite?

This script answers one question only: which (calculator, pbc) combinations
accept the 96 atom canonical slab and return an energy and forces. It performs
NO relaxation and writes NO structures. It is safe and free to run.

Usage:
    python scripts/pbc_feasibility_probe.py \
        --slab inputs/slab_clean_2x2_canonical.in \
        --out outputs/pbc_feasibility_probe.json
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
import traceback

import numpy as np
from ase.io import read


def package_versions() -> dict:
    """Record every version that could change the numbers below."""
    vers = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
    }
    for name, mod in (("numpy", "numpy"), ("ase", "ase"), ("xtb", "xtb"),
                      ("tblite", "tblite")):
        try:
            vers[name] = __import__(mod).__version__
        except Exception as exc:  # noqa: BLE001
            vers[name] = f"UNAVAILABLE: {type(exc).__name__}: {exc}"
    # The compiled xtb library version is distinct from the xtb-python wrapper
    # version; both matter for reproducibility.
    try:
        from xtb.libxtb import VERSION as XTB_LIB_VERSION
        vers["xtb_library"] = str(XTB_LIB_VERSION)
    except Exception as exc:  # noqa: BLE001
        vers["xtb_library"] = f"UNAVAILABLE: {type(exc).__name__}: {exc}"
    return vers


def exception_chain(exc: BaseException) -> list:
    """Every message in a chained exception, outermost first.

    ASE wraps a calculator failure in its own CalculationFailed, so str(exc)
    returns "xtb could not evaluate input" and the library's own diagnosis,
    which is the entire point of this probe, sits in __context__ one level
    down. Walking the chain is what surfaces "Multipoles not available with
    PBC" instead of a generic wrapper message.
    """
    msgs, seen, cur = [], set(), exc
    while cur is not None and id(cur) not in seen:
        seen.add(id(cur))
        msgs.append(f"{type(cur).__name__}: {cur}".strip())
        cur = cur.__cause__ or cur.__context__
    return msgs


def make_calculator(kind: str, method: str = "GFN2-xTB"):
    """Return a fresh calculator instance. Fresh per trial, so that a failure
    in one trial cannot leave cached state behind in the next."""
    if kind == "xtb-python":
        from xtb.ase.calculator import XTB
        return XTB(method=method)
    if kind == "tblite":
        from tblite.ase import TBLite
        return TBLite(method=method)
    raise ValueError(f"unknown calculator kind: {kind}")


def describe_structure(atoms) -> dict:
    z = atoms.positions[:, 2]
    cell = atoms.get_cell()
    symbols = atoms.get_chemical_symbols()
    return {
        "n_atoms": len(atoms),
        "formula": atoms.get_chemical_formula(),
        "n_Ir": symbols.count("Ir"),
        "n_O": symbols.count("O"),
        "n_H": symbols.count("H"),
        "cell_lengths_A": [float(v) for v in cell.lengths()],
        "cell_angles_deg": [float(v) for v in cell.angles()],
        "z_min_A": float(z.min()),
        "z_max_A": float(z.max()),
        "slab_thickness_A": float(z.max() - z.min()),
        "vacuum_gap_along_c_A": float(cell.lengths()[2] - (z.max() - z.min())),
    }


def trial(atoms_template, kind: str, pbc, method: str) -> dict:
    """Attempt one single point. Never raises: a failure is a result."""
    atoms = atoms_template.copy()
    atoms.set_pbc(pbc)
    rec = {
        "calculator": kind,
        "method": method,
        "pbc": [bool(p) for p in atoms.get_pbc()],
    }
    t0 = time.perf_counter()
    try:
        atoms.calc = make_calculator(kind, method)
        energy = float(atoms.get_potential_energy())
        forces = np.asarray(atoms.get_forces(), dtype=float)
        fnorm = np.linalg.norm(forces, axis=1)
        rec.update({
            "status": "ok",
            "energy_eV": energy,
            "fmax_all_atoms_eV_per_A": float(fnorm.max()),
            "fmax_argmax_index": int(fnorm.argmax()),
            "frms_all_atoms_eV_per_A": float(np.sqrt((forces ** 2).sum(axis=1).mean())),
            "force_sum_eV_per_A": [float(v) for v in forces.sum(axis=0)],
        })
    except BaseException as exc:  # noqa: BLE001  (xtb can raise non-Exception)
        chain = exception_chain(exc)
        rec.update({
            "status": "failed",
            "error_type": type(exc).__name__,
            "error": str(exc),
            # The library's own message, which the ASE wrapper hides.
            "error_chain": chain,
            "error_root": chain[-1] if chain else str(exc),
            "traceback_tail": traceback.format_exc().strip().splitlines()[-4:],
        })
    rec["wall_time_s"] = round(time.perf_counter() - t0, 3)
    return rec


PBC_CHOICES = {
    "FFF": (False, False, False),   # current pipeline: finite cluster
    "TTF": (True, True, False),     # what the request asks for: 2D slab
    "TTT": (True, True, True),      # 3D, relying on the 42.013 A c axis
}


def run_trial_subprocess(slab: str, kind: str, pbc_key: str, method: str) -> dict:
    """Run one trial in a separate interpreter.

    A tight-binding library that aborts on an unsupported boundary condition
    takes the whole process down with it, so an in-process loop would lose
    every later trial. Isolating each trial turns a crash into a recorded
    result: returncode -11 is SIGSEGV, -6 is SIGABRT.
    """
    import subprocess
    proc = subprocess.run(
        [sys.executable, __file__, "--single", "--slab", slab,
         "--kind", kind, "--pbc-key", pbc_key, "--method", method],
        capture_output=True, text=True, timeout=3600,
    )
    for line in proc.stdout.splitlines():
        if line.startswith("@@RESULT@@"):
            return json.loads(line[len("@@RESULT@@"):])
    # No result line means the child died before it could report.
    rc = proc.returncode
    signal_name = ""
    if rc < 0:
        import signal as _sig
        try:
            signal_name = _sig.Signals(-rc).name
        except ValueError:
            signal_name = f"signal {-rc}"
    return {
        "calculator": kind,
        "method": method,
        "pbc": [bool(p) for p in PBC_CHOICES[pbc_key]],
        "status": "crashed",
        "returncode": rc,
        "signal": signal_name,
        "stderr_tail": proc.stderr.strip().splitlines()[-6:],
        "stdout_tail": proc.stdout.strip().splitlines()[-6:],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--slab", default="inputs/slab_clean_2x2_canonical.in")
    ap.add_argument("--out", default="outputs/pbc_feasibility_probe.json")
    ap.add_argument("--method", default="GFN2-xTB")
    ap.add_argument("--calculators", nargs="+", default=["xtb-python", "tblite"],
                choices=["xtb-python", "tblite"],
                help="Engines to trial. CI passes --calculators tblite, because "
                    "tblite-python 0.5.0 cannot share an environment with xtb.")
    ap.add_argument("--single", action="store_true",
                    help="Internal: run one trial and print a @@RESULT@@ line.")
    ap.add_argument("--kind", default="xtb-python")
    ap.add_argument("--pbc-key", default="FFF", choices=sorted(PBC_CHOICES))
    args = ap.parse_args()

    if args.single:
        atoms0 = read(args.slab, format="espresso-in")
        rec = trial(atoms0, args.kind, PBC_CHOICES[args.pbc_key], args.method)
        print("@@RESULT@@" + json.dumps(rec), flush=True)
        return

    # ASE reads a Quantum ESPRESSO input as fully 3D periodic. Every trial sets
    # its own pbc explicitly, so the value here is only the starting point.
    atoms0 = read(args.slab, format="espresso-in")

    report = {
        "probe": "step1_pbc_feasibility",
        "slab_file": args.slab,
        "method": args.method,
        "versions": package_versions(),
        "structure": describe_structure(atoms0),
        "pbc_as_read_from_qe_input": [bool(p) for p in atoms0.get_pbc()],
        "trials": [],
    }

    # GFN1-xTB is included as a control, not as a proposed production method.
    # It isolates "this Hamiltonian cannot do PBC" from "this build is broken":
    # if GFN1 runs periodically where GFN2 does not, the obstacle is the GFN2
    # multipole electrostatics, not the periodic machinery or the install.
    for kind, meth in (("xtb-python", "GFN2-xTB"), ("tblite", "GFN2-xTB"),
                       ("xtb-python", "GFN1-xTB"), ("tblite", "GFN1-xTB")):
        if kind not in args.calculators:
            continue
        for pbc_key in ("FFF", "TTF", "TTT"):
            rec = run_trial_subprocess(args.slab, kind, pbc_key, meth)
            report["trials"].append(rec)
            detail = rec.get("energy_eV")
            if detail is None:
                detail = rec.get("error_type") or rec.get("signal") or ""
            print(f"[{rec['status']:>7}] {kind:<11} {meth} pbc={pbc_key}  {detail}",
                  flush=True)

    from pathlib import Path
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(f"[ok] wrote {out}")


if __name__ == "__main__":
    main()
