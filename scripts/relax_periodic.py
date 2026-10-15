#!/usr/bin/env python3
"""Steps 3 and 4: relaxation under a chosen boundary condition, with the
frozen optimiser configuration otherwise unchanged.

EXECUTION STATUS: NOT RUN. This script has been written and reviewed but has
not been executed. Every number it prints will be produced on your machine.

What is held fixed, read from inputs/frozen_optimiser_config.json and not
overridden here: GFN2-xTB by default, ASE BFGS, fmax 0.05 eV/A, 200 steps,
FixAtoms on z below z_freeze_A. The single deliberate change is the boundary
condition, supplied by --pbc.

Reports for every run:
  - convergence and step count
  - minimum O-O distance, minimum Ir-O distance
  - whether any O-O contact below 1.60 A forms
  - maximum displacement of mobile atoms

Additionally, when --site is given (step 4):
  - H to site O distance
  - identity of the nearest O to H at every frame of the trajectory
  - desorption check: nearest O distance above 2.5 A and increasing

All distances use the minimum image convention whenever the corresponding
pbc flag is set, so the periodic and open runs are measured consistently.

This script never writes to site_intended and never generates a new site.

Usage, step 3 (clean slab):
    python scripts/relax_periodic.py --structure outputs/clean/slab_clean_ready.traj \
        --pbc TTT --out outputs/step3_clean_periodic

Usage, step 4 (H* on o69, only if step 3 is clean):
    python scripts/relax_periodic.py --structure outputs/H_star_o69/..._ready.traj \
        --pbc TTT --site 69 --out outputs/step4_h_o69_periodic
"""
from __future__ import annotations

import argparse
import json
import traceback
from pathlib import Path

import numpy as np
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.optimize import BFGS

PBC_CHOICES = {
    "FFF": (False, False, False),
    "TTF": (True, True, False),
    "TTT": (True, True, True),
}


def exception_chain(exc: BaseException) -> list:
    """Every message in a chained exception, outermost first. ASE hides the
    calculator library's own diagnosis one level down in __context__, so
    str(exc) alone returns a generic wrapper message. A periodic GFN2
    attempt fails on "Multipoles not available with PBC", and that string
    must survive into the results JSON."""
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


def pair_minima(atoms) -> dict:
    """Minimum O-O and Ir-O separations, honouring periodicity where set."""
    d = atoms.get_all_distances(mic=any(atoms.get_pbc()))
    syms = np.array(atoms.get_chemical_symbols())
    iu = np.triu_indices(len(atoms), k=1)
    pair = [tuple(sorted((syms[i], syms[j]))) for i, j in zip(*iu)]
    pair = np.array(["-".join(p) for p in pair])
    dv = d[iu]
    out = {}
    for label, key in (("min_OO_A", "O-O"), ("min_IrO_A", "Ir-O"),
                       ("min_IrIr_A", "Ir-Ir")):
        sel = pair == key
        if sel.any():
            k = int(np.argmin(dv[sel]))
            idx_pairs = np.array(list(zip(*iu)))[sel]
            out[label] = float(dv[sel][k])
            out[label.replace("_A", "_pair")] = [int(v) for v in idx_pairs[k]]
        else:
            out[label] = None
    return out


def short_OO_contacts(atoms, cutoff: float = 1.60) -> list:
    d = atoms.get_all_distances(mic=any(atoms.get_pbc()))
    syms = np.array(atoms.get_chemical_symbols())
    hits = []
    oi = np.where(syms == "O")[0]
    for a_i in range(len(oi)):
        for b_i in range(a_i + 1, len(oi)):
            i, j = int(oi[a_i]), int(oi[b_i])
            if d[i, j] < cutoff:
                hits.append({"i": i, "j": j, "distance_A": float(d[i, j])})
    return sorted(hits, key=lambda h: h["distance_A"])


def hydrogen_tracking(traj_path: Path, site_index: int) -> dict:
    """Nearest O to H, frame by frame, plus a desorption verdict."""
    frames = list(Trajectory(str(traj_path)))
    if not frames:
        return {"error": "empty trajectory"}
    syms = np.array(frames[0].get_chemical_symbols())
    h_idx = np.where(syms == "H")[0]
    if len(h_idx) != 1:
        return {"error": f"expected exactly one H, found {len(h_idx)}"}
    h = int(h_idx[0])
    o_idx = np.where(syms == "O")[0]

    series = []
    for n, fr in enumerate(frames):
        mic = any(fr.get_pbc())
        d_h = fr.get_distances(h, o_idx, mic=mic)
        k = int(np.argmin(d_h))
        series.append({
            "frame": n,
            "nearest_O_index": int(o_idx[k]),
            "nearest_O_distance_A": float(d_h[k]),
            "d_H_to_site_O_A": float(fr.get_distance(h, site_index, mic=mic)),
        })

    last = series[-1]
    # Desorbing means both far and still receding. A single far frame is not
    # enough, so the trend over the closing frames is checked as well.
    tail = [s["nearest_O_distance_A"] for s in series[-min(10, len(series)):]]
    increasing = len(tail) > 1 and (tail[-1] - tail[0]) > 0.0
    return {
        "H_index": h,
        "site_O_index": site_index,
        "final_d_H_to_site_O_A": last["d_H_to_site_O_A"],
        "final_nearest_O_index": last["nearest_O_index"],
        "final_nearest_O_distance_A": last["nearest_O_distance_A"],
        "site_O_is_nearest_at_end": last["nearest_O_index"] == site_index,
        "nearest_O_changed_during_run": len(
            {s["nearest_O_index"] for s in series}) > 1,
        "nearest_O_identities_seen": sorted(
            {s["nearest_O_index"] for s in series}),
        "desorbing": bool(last["nearest_O_distance_A"] > 2.5 and increasing),
        "series": series,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--structure", required=True,
                    help="Input structure (.traj or any ASE readable format)")
    ap.add_argument("--config", default="inputs/frozen_optimiser_config.json")
    ap.add_argument("--pbc", default="TTT", choices=sorted(PBC_CHOICES))
    ap.add_argument("--calculator", default="xtb-python",
                    choices=["xtb-python", "tblite"])
    ap.add_argument("--method", default=None,
                    help="Overrides gfn_level from the frozen config")
    ap.add_argument("--site", type=int, default=None,
                    help="Zero based ASE index of the site O, step 4 only")
    ap.add_argument("--out", required=True, help="Output directory")
    args = ap.parse_args()

    cfg = json.loads(Path(args.config).read_text())
    method = args.method or cfg["gfn_level"]
    fmax = float(cfg["fmax_target_eV_per_A"])
    steps = int(cfg["max_steps"])
    z_freeze = float(cfg["constraint"]["z_freeze_A"])

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    base = Path(args.structure).stem
    results_path = outdir / f"{base}_{args.pbc}_results.json"

    result = {
        "structure_file": args.structure,
        "calculator": args.calculator,
        "method": method,
        "pbc": args.pbc,
        "fmax_target": fmax,
        "max_steps": steps,
        "z_freeze_A": z_freeze,
        "config_name": cfg.get("config_name"),
        "status": "started",
    }
    results_path.write_text(json.dumps(result, indent=2))

    try:
        atoms = read(args.structure)
        atoms.set_pbc(PBC_CHOICES[args.pbc])
        mask = atoms.positions[:, 2] < z_freeze
        atoms.set_constraint(FixAtoms(mask=mask))
        pos0 = atoms.positions.copy()

        atoms.calc = make_calculator(args.calculator, method)
        result.update({
            "n_atoms": len(atoms),
            "n_frozen": int(mask.sum()),
            "n_free": int((~mask).sum()),
            "initial": {**pair_minima(atoms),
                        "energy_eV": float(atoms.get_potential_energy())},
        })

        logfile = outdir / f"{base}_{args.pbc}_opt.log"
        trajlog = outdir / f"{base}_{args.pbc}_opt.traj"
        opt = BFGS(atoms, logfile=str(logfile), trajectory=str(trajlog))
        print(f"[opt] {base} pbc={args.pbc}: start", flush=True)
        opt.run(fmax=fmax, steps=steps)
        print(f"[opt] {base} pbc={args.pbc}: end", flush=True)

        f = np.asarray(atoms.get_forces(apply_constraint=False))
        fn = np.linalg.norm(f, axis=1)
        disp = np.linalg.norm(atoms.positions - pos0, axis=1)
        contacts = short_OO_contacts(atoms, 1.60)

        result.update({
            "status": "ok",
            "n_steps": int(getattr(opt, "nsteps", 0)),
            "converged": bool(fn[~mask].max() <= fmax),
            "final": {
                **pair_minima(atoms),
                "energy_eV": float(atoms.get_potential_energy()),
                "fmax_free_atoms_eV_per_A": float(fn[~mask].max()),
            },
            "max_displacement_mobile_A": float(disp[~mask].max()),
            "max_displacement_mobile_index": int(np.argmax(np.where(~mask, disp, -1))),
            "short_OO_below_1.60A_count": len(contacts),
            "short_OO_below_1.60A": contacts[:20],
            "final_structure": str(outdir / f"{base}_{args.pbc}_final.traj"),
        })
        write(outdir / f"{base}_{args.pbc}_final.traj", atoms)

        if args.site is not None:
            result["hydrogen"] = hydrogen_tracking(trajlog, args.site)

    except BaseException as exc:  # noqa: BLE001
        chain = exception_chain(exc)
        result.update({
            "status": "error",
            "error_type": type(exc).__name__,
            "error": str(exc),
            "error_chain": chain,
            "error_root": chain[-1] if chain else str(exc),
            "traceback": traceback.format_exc(),
        })
    finally:
        results_path.write_text(json.dumps(result, indent=2))
        print(f"[ok] wrote {results_path}", flush=True)


if __name__ == "__main__":
    main()
