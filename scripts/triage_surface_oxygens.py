#!/usr/bin/env python3
"""Triage the surface oxygens of an IrO2 slab before any state generation.

This is the second gate in the pipeline (the first is the bond-length
geometry audit, which this script also runs). It answers three questions
and refuses to bless a slab unless all three come back right:

1. Is the local geometry rutile? (bond-length audit: Ir-O first shell in
   1.90 to 2.10 A, no O-O contact below 2.2 A)
2. Which oxygens are surface oxygens, and what kind of site is each one?
   (terminal 1-fold, bridging 2-fold, in-plane 3-fold, by Ir coordination)
3. Does the surface match the canonical stoichiometric rutile (110)
   termination? (bridging 2-fold O rows above a plane of 5-fold cus and
   6-fold Ir with 3-fold in-plane O; no 1-fold terminal O anywhere)

Site indices are coordinates in a file, not identities: after any change
to the slab, targets must be re-chosen from this script's output by site
type, never carried over by index.

Usage:
    python scripts/triage_surface_oxygens.py SLAB [--json out.json]
        [--surface-threshold 1.5] [--ir-cutoff 2.2] [--no-verdict]

SLAB may be a .traj file or a Quantum ESPRESSO input (.in).
Exit code 0 only if the audit passes and (unless --no-verdict) the
termination is canonical (110). Workflows can gate on the exit code.

Requires: ase, numpy. No other dependencies.
"""

import argparse
import json
import sys

import numpy as np
from ase.io import read

SITE_NAMES = {
    1: "terminal (1-fold)",
    2: "bridging (2-fold)",
    3: "in-plane (3-fold)",
}


def load_slab(path):
    if str(path).endswith(".in"):
        return read(path, format="espresso-in")
    return read(path)


def bond_audit(atoms):
    # In-plane periodicity must be respected when counting neighbours,
    # otherwise atoms at the cell edge lose bonds across the boundary and
    # coordination-based checks misfire. The optimisation stage stores
    # slabs with pbc switched off (xTB convention), so restore x/y
    # periodicity from the cell for analysis only.
    atoms = atoms.copy()
    atoms.set_pbc((True, True, False))
    sym = np.array(atoms.get_chemical_symbols())
    ir = np.where(sym == "Ir")[0]
    ox = np.where(sym == "O")[0]
    D = atoms.get_all_distances(mic=True)
    min_iro = float(min(D[i][ox].min() for i in ir))
    oo = D[np.ix_(ox, ox)]
    iu = np.triu_indices_from(oo, k=1)
    min_oo = float(oo[iu].min())
    n_close = int((oo[iu] < 2.2).sum())
    passed = 1.90 <= min_iro <= 2.10 and n_close == 0
    return {
        "min_Ir_O_A": round(min_iro, 3),
        "min_O_O_A": round(min_oo, 3),
        "O_O_contacts_below_2p2_A": n_close,
        "passed": bool(passed),
    }, D, ir, ox


def triage(atoms, surface_threshold, ir_cutoff):
    audit, D, ir, ox = bond_audit(atoms)
    z = atoms.get_positions()[:, 2]
    z_top = float(z.max())

    surface_o = []
    for i in ox:
        if (z_top - z[i]) < surface_threshold:
            n_ir = int((D[i][ir] < ir_cutoff).sum())
            surface_o.append({
                "index": int(i),
                "label": f"o{int(i)}",
                "z_A": round(float(z[i]), 3),
                "depth_below_top_A": round(z_top - float(z[i]), 3),
                "ir_coordination": n_ir,
                "site_type": SITE_NAMES.get(n_ir, f"{n_ir}-fold"),
            })

    # Top-layer Ir plane: the Ir closest to the surface oxygens.
    ir_z = z[ir]
    top_ir_z = float(ir_z.max())
    top_ir = [i for i in ir if abs(z[i] - top_ir_z) < 0.1]
    top_ir_coord = sorted(int((D[i][ox] < ir_cutoff).sum()) for i in top_ir)

    # Canonical stoichiometric rutile (110) verdict.
    types = [s["ir_coordination"] for s in surface_o]
    reasons = []
    if any(t == 1 for t in types):
        n1 = sum(1 for t in types if t == 1)
        reasons.append(f"{n1} terminal 1-fold O in the surface region "
                       "(canonical (110) has none)")
    if not any(t == 2 for t in types):
        reasons.append("no 2-fold bridging O found "
                       "(canonical (110) is terminated by bridging O rows)")
    if not set(top_ir_coord) <= {5, 6}:
        reasons.append(f"top-plane Ir coordinations {sorted(set(top_ir_coord))} "
                       "(canonical (110) exposes only 5-fold cus and 6-fold Ir)")
    canonical = not reasons

    return {
        "n_atoms": len(atoms),
        "z_top_A": round(z_top, 3),
        "surface_threshold_A": surface_threshold,
        "ir_cutoff_A": ir_cutoff,
        "bond_audit": audit,
        "surface_oxygens": surface_o,
        "bridging_o_indices": [s["index"] for s in surface_o
                               if s["ir_coordination"] == 2],
        "in_plane_o_indices": [s["index"] for s in surface_o
                               if s["ir_coordination"] == 3],
        "terminal_o_indices": [s["index"] for s in surface_o
                               if s["ir_coordination"] == 1],
        "top_ir_plane_coordinations": top_ir_coord,
        "canonical_110": bool(canonical),
        "canonical_110_failures": reasons,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("slab", help="slab structure (.traj or espresso .in)")
    p.add_argument("--surface-threshold", type=float, default=1.5,
                   help="O within this depth of z_top counts as surface "
                        "(same rule as iro2_slab_setup.py)")
    p.add_argument("--ir-cutoff", type=float, default=2.2,
                   help="Ir neighbour cutoff for coordination counting")
    p.add_argument("--json", help="also write the full report to this path")
    p.add_argument("--no-verdict", action="store_true",
                   help="report only; do not fail on a non-canonical "
                        "termination (bond audit still gates)")
    args = p.parse_args()

    atoms = load_slab(args.slab)
    report = triage(atoms, args.surface_threshold, args.ir_cutoff)

    a = report["bond_audit"]
    print(f"slab: {args.slab} ({report['n_atoms']} atoms)")
    print(f"bond audit: min Ir-O {a['min_Ir_O_A']} A, min O-O {a['min_O_O_A']} A, "
          f"O-O<2.2A {a['O_O_contacts_below_2p2_A']} "
          f"-> {'PASS' if a['passed'] else 'FAIL'}")
    print(f"\nsurface oxygens (within {args.surface_threshold} A of "
          f"z_top={report['z_top_A']}):")
    print(f"{'label':6s} {'z (A)':>8s} {'depth':>7s} {'Ir-coord':>9s}  site type")
    for s in report["surface_oxygens"]:
        print(f"{s['label']:6s} {s['z_A']:8.2f} {s['depth_below_top_A']:7.2f} "
              f"{s['ir_coordination']:9d}  {s['site_type']}")
    print(f"\nbridging O (H* candidate sites): "
          f"{[f'o{i}' for i in report['bridging_o_indices']] or 'NONE'}")
    print(f"in-plane O: {[f'o{i}' for i in report['in_plane_o_indices']] or 'none'}")
    if report["terminal_o_indices"]:
        print(f"terminal O (should not exist on canonical (110)): "
              f"{[f'o{i}' for i in report['terminal_o_indices']]}")
    print(f"top Ir plane coordinations: {report['top_ir_plane_coordinations']}")

    if report["canonical_110"]:
        print("\ntermination: CANONICAL rutile (110)")
    else:
        print("\ntermination: NOT canonical rutile (110):")
        for r in report["canonical_110_failures"]:
            print(f"  - {r}")

    if args.json:
        with open(args.json, "w") as f:
            json.dump(report, f, indent=2)
        print(f"\nreport written to {args.json}")

    if not a["passed"]:
        sys.exit("FAIL: bond-length audit")
    if not report["canonical_110"] and not args.no_verdict:
        sys.exit("FAIL: non-canonical termination; do not generate states "
                 "from this slab (use --no-verdict to report only)")
    print("\nGATE PASSED: states may be generated from this slab.")


if __name__ == "__main__":
    main()
