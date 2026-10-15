#!/usr/bin/env python3
"""Is pbc (True, True, False) honoured as genuine 2D periodicity, or is it
silently promoted to 3D using the c lattice vector?

The feasibility probe found that (T,T,F) and (T,T,T) give identical energies.
That is consistent with two different explanations:

  (a) the calculator implements true 2D periodicity, and the 30 A vacuum makes
      the 3D answer numerically indistinguishable from it, or
  (b) the calculator ignores the third pbc flag and always uses the full 3D
      lattice, so (T,T,F) is 3D in disguise.

These are distinguished by changing c alone. Under (a) the (T,T,F) energy is
independent of c. Under (b) it moves with c. Nothing else in the structure
changes: the atoms are untouched and only the third lattice vector is scaled.

Usage:
    python scripts/pbc_2d_vs_3d_check.py --slab inputs/slab_clean_2x2_canonical.in
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ase.io import read


def exception_chain(exc: BaseException) -> list:
    """Every message in a chained exception, outermost first. ASE hides the
    calculator library's own diagnosis one level down in __context__."""
    msgs, seen, cur = [], set(), exc
    while cur is not None and id(cur) not in seen:
        seen.add(id(cur))
        msgs.append(f"{type(cur).__name__}: {cur}".strip())
        cur = cur.__cause__ or cur.__context__
    return msgs


def energy_at_c(slab: str, kind: str, method: str, pbc, c_length: float) -> dict:
    rec = {
        "calculator": kind,
        "method": method,
        "pbc": [bool(p) for p in pbc],
        "c_A": float(c_length),
    }
    try:
        atoms = read(slab, format="espresso-in")
        cell = atoms.get_cell().array.copy()
        cell[2, 2] = c_length
        # scale_atoms=False: atomic positions are held fixed in Cartesian
        # space, so the only thing that changes between trials is the vacuum.
        atoms.set_cell(cell, scale_atoms=False)
        atoms.set_pbc(pbc)
        if kind == "xtb-python":
            from xtb.ase.calculator import XTB
            atoms.calc = XTB(method=method)
        else:
            from tblite.ase import TBLite
            atoms.calc = TBLite(method=method)
        rec["energy_eV"] = float(atoms.get_potential_energy())
        rec["status"] = "ok"
    except BaseException as exc:  # noqa: BLE001
        # A refusal is a result. With --method GFN2-xTB the xtb-python rows
        # are expected to fail on "Multipoles not available with PBC", and
        # that must not stop the tblite rows from being measured.
        chain = exception_chain(exc)
        rec["status"] = "failed"
        rec["energy_eV"] = None
        rec["error_type"] = type(exc).__name__
        rec["error"] = str(exc)
        rec["error_chain"] = chain
        rec["error_root"] = chain[-1] if chain else str(exc)
    return rec


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--slab", default="inputs/slab_clean_2x2_canonical.in")
    ap.add_argument("--out", default="outputs/pbc_2d_vs_3d_check.json")
    # GFN1-xTB is the only method that runs periodically at all, so it is the
    # only one this question can be asked of.
    ap.add_argument("--method", default="GFN1-xTB")
    ap.add_argument("--calculators", nargs="+", default=["xtb-python", "tblite"],
                    choices=["xtb-python", "tblite"],
                    help="Restrict to one engine, e.g. --calculators tblite")
    ap.add_argument("--c0", type=float, default=42.01316398881907,
                    help="Canonical c, giving 30 A vacuum")
    ap.add_argument("--c1", type=float, default=62.01316398881907,
                    help="Comparison c, atoms unmoved")
    args = ap.parse_args()

    c0, c1 = args.c0, args.c1

    rows = []
    for kind in args.calculators:
        for pbc in ((True, True, False), (True, True, True)):
            for c in (c0, c1):
                rec = energy_at_c(args.slab, kind, args.method, pbc, c)
                rows.append(rec)
                lbl = "".join("T" if p else "F" for p in pbc)
                if rec["status"] == "ok":
                    print(f"[ok]     {kind:<11} pbc={lbl} c={c:7.3f}  "
                          f"E={rec['energy_eV']:.9f} eV", flush=True)
                else:
                    print(f"[failed] {kind:<11} pbc={lbl} c={c:7.3f}  "
                          f"{rec['error_root']}", flush=True)

    # Summarise the discriminator: dE under a change in vacuum. A pair with
    # any failed member cannot yield a difference, and is reported as such.
    summary = []
    for kind in args.calculators:
        for pbc_lbl in ("TTF", "TTT"):
            sel = [r for r in rows
                   if r["calculator"] == kind
                   and "".join("T" if p else "F" for p in r["pbc"]) == pbc_lbl]
            lo = next(r for r in sel if abs(r["c_A"] - c0) < 1e-9)
            hi = next(r for r in sel if abs(r["c_A"] - c1) < 1e-9)
            if lo["status"] != "ok" or hi["status"] != "ok":
                summary.append({"calculator": kind, "pbc": pbc_lbl,
                                "dE_for_vacuum_change_eV": None,
                                "verdict": "not determined: a single point failed"})
                print(f"  -> {kind:<11} {pbc_lbl}: not determined, "
                      f"a single point failed")
                continue
            d = hi["energy_eV"] - lo["energy_eV"]
            summary.append({
                "calculator": kind, "pbc": pbc_lbl,
                "dE_for_vacuum_change_eV": d,
                "verdict": ("c is used: this is 3D periodicity"
                            # 1e-4 eV, not 1e-6. SCF convergence noise on a
                            # 96-atom slab is of order 1e-6 eV, so a 1e-6
                            # threshold let tblite's flat response flip
                            # between verdicts on noise alone (Linux CI gave
                            # -2.3e-7 for TTF and +1.0e-6 for TTT). xtb's
                            # 0.718 eV response is four orders above this.
                            if abs(d) > 1e-4 else
                            "c is irrelevant: consistent with true 2D"),
            })
            print(f"  -> {kind:<11} {pbc_lbl}: dE = {d:+.9f} eV  "
                  f"({summary[-1]['verdict']})")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"c_values_A": [c0, c1],
                               "trials": rows, "summary": summary}, indent=2))
    print(f"[ok] wrote {out}")


if __name__ == "__main__":
    main()
