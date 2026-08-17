#!/usr/bin/env python3
"""Build a canonical stoichiometric rutile IrO2(110) 2x2 slab, gate-verified.

Construction: the bulk is built from published experimental rutile IrO2
parameters (Bolzan et al. 1997: a=4.4983 A, c=3.1544 A, u=0.3058, space
group P4_2/mnm). A z-periodic (110) stack is generated, rolled through its
periodic boundary in fine steps, and each candidate cut is judged by the
surface-oxygen triage gate (scripts/triage_surface_oxygens.py). Only cuts
that pass the bond audit AND the canonical-(110) verdict survive; the
middle of the widest passing window is selected for robustness. The chosen
slab is gated one final time before anything is written.

ASE's default surface() cut lands mid-trilayer and produces an oxygen-rich
top face (terminal O capping the cus Ir rows); this script exists because
of that. See the repository README for the correction history.

Usage:
    python scripts/build_canonical_110_slab.py --out inputs/slab_clean_2x2_canonical
Writes <out>.traj and <out>.in and prints the surface-site table and the
recommended z_freeze for the frozen optimiser config.

Requires: ase, numpy, and scripts/triage_surface_oxygens.py.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from ase.build import surface
from ase.io import write
from ase.spacegroup import crystal

sys.path.insert(0, str(Path(__file__).resolve().parent))
from triage_surface_oxygens import triage  # noqa: E402

A, C, U = 4.4983, 3.1544, 0.3058  # Bolzan et al. 1997, rutile IrO2


def build_candidate(periodic, roll, repeat, vacuum):
    s = periodic.copy()
    sp = s.get_scaled_positions()
    sp[:, 2] = (sp[:, 2] + roll) % 1.0
    s.set_scaled_positions(sp)
    s = s.repeat((repeat[0], repeat[1], 1))
    s.center(vacuum=vacuum, axis=2)
    s.set_pbc((True, True, False))
    return s


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--layers", type=int, default=4)
    p.add_argument("--repeat", type=int, nargs=2, default=[2, 2])
    p.add_argument("--vacuum", type=float, default=15.0)
    p.add_argument("--scan-steps", type=int, default=96)
    p.add_argument("--out", default="inputs/slab_clean_2x2_canonical")
    args = p.parse_args()

    bulk = crystal(["Ir", "O"], basis=[(0, 0, 0), (U, U, 0)],
                   spacegroup=136, cellpar=[A, A, C, 90, 90, 90])
    periodic = surface(bulk, (1, 1, 0), layers=args.layers,
                       vacuum=None, periodic=True)

    rolls = np.linspace(0, 1, args.scan_steps, endpoint=False)
    passing = []
    for d in rolls:
        r = triage(build_candidate(periodic, d, args.repeat, args.vacuum),
                   1.5, 2.2)
        if r["bond_audit"]["passed"] and r["canonical_110"]:
            passing.append(d)
    if not passing:
        sys.exit("No roll produced a canonical (110) termination; "
                 "inspect the stack before proceeding.")

    # Middle of the widest contiguous passing window (wrap-aware).
    step = 1.0 / args.scan_steps
    windows, cur = [], [passing[0]]
    for d in passing[1:]:
        if abs(d - cur[-1]) <= step * 1.5:
            cur.append(d)
        else:
            windows.append(cur)
            cur = [d]
    windows.append(cur)
    if len(windows) > 1 and abs((passing[0] + 1.0) - passing[-1]) <= step * 1.5:
        windows[0] = windows.pop() + windows[0]  # wrap-around join
    best = max(windows, key=len)
    roll = best[len(best) // 2]

    slab = build_candidate(periodic, roll, args.repeat, args.vacuum)
    report = triage(slab, 1.5, 2.2)
    assert report["bond_audit"]["passed"] and report["canonical_110"], \
        "final gate failed; refusing to write"

    z = slab.get_positions()[:, 2]
    z_freeze = (z.max() + z.min()) / 2

    print(f"passing rolls: {len(passing)}/{args.scan_steps}; "
          f"selected roll {roll:.4f} (window of {len(best)})")
    print(f"atoms: {len(slab)}")
    print(f"bridging O (H* candidate sites): "
          f"{[f'o{i}' for i in report['bridging_o_indices']]}")
    print(f"in-plane O: {[f'o{i}' for i in report['in_plane_o_indices']]}")
    print(f"top Ir plane coordinations: "
          f"{sorted(set(report['top_ir_plane_coordinations']))}")
    print(f"recommended z_freeze: {z_freeze:.2f} A "
          "(update the frozen optimiser config, new version, new hash)")

    write(f"{args.out}.traj", slab)
    write(f"{args.out}.in", slab, format="espresso-in",
          pseudopotentials={"Ir": "Ir.us.z_31.ld1.psl.v1.0.0-high.upf",
                            "O": "O.pbe-n-kjpaw_psl.0.1.UPF"})
    print(f"wrote {args.out}.traj and {args.out}.in (gate-verified)")


if __name__ == "__main__":
    main()
