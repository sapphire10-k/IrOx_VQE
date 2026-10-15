#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import traceback
from pathlib import Path

from ase.io import read, write
from ase.optimize import BFGS
from xtb.ase.calculator import XTB


def infer_base_name(traj_path: str) -> str:
    return Path(traj_path).stem


def compute_fmax(atoms) -> float:
    forces = atoms.get_forces()
    # max norm of force vectors
    return float(((forces**2).sum(axis=1) ** 0.5).max())


def main() -> None:
    ap = argparse.ArgumentParser(description="Optimize an ASE .traj structure with xTB and write robust results JSON.")
    ap.add_argument("traj", help="Input .traj file")
    ap.add_argument("--output", required=True, help="Output directory")
    ap.add_argument("--fmax", type=float, default=0.05, help="Convergence threshold (eV/Å)")
    ap.add_argument("--steps", type=int, default=200, help="Maximum optimizer steps")
    ap.add_argument("--method", default="GFN2-xTB", help="xTB method")
    ap.add_argument("--config", default="inputs/frozen_optimiser_config.json",
                    help="Frozen optimiser config, checked against the structure")
    args = ap.parse_args()

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    base_name = infer_base_name(args.traj)
    results_path = outdir / f"{base_name}_results.json"

    # Write a "header" immediately so we always have a file even if something crashes early.
    result = {
        "base_name": base_name,
        "structure_file": str(Path(args.traj)),
        "method": args.method,
        "fmax_target": float(args.fmax),
        "max_steps": int(args.steps),
        "status": "started",
    }
    results_path.write_text(json.dumps(result, indent=2))

    try:
        atoms = read(args.traj)
        cfg = json.loads(Path(args.config).read_text())
        expected = tuple(cfg["pbc"])
        got = tuple(bool(p) for p in atoms.get_pbc())
        if got != expected:
            raise SystemExit(
                f"{args.traj}: pbc {got} does not match frozen config "
                f"{cfg['config_name']} which requires {expected}. "
                f"Rebuild the states with the current generate_states.py.")
        atoms.calc = XTB(method=args.method)

        # initial diagnostics
        e_init = float(atoms.get_potential_energy())
        f_init = compute_fmax(atoms)
        result.update({"e_init": e_init, "fmax_init": f_init})

        logfile = outdir / f"{base_name}_opt.log"
        trajlog = outdir / f"{base_name}_opt.traj"
        opt = BFGS(atoms, logfile=str(logfile), trajectory=str(trajlog))
        print(f"[opt] {base_name}: start (fmax={args.fmax}, steps={args.steps})", flush=True)
        opt.run(fmax=args.fmax, steps=args.steps)
        print(f"[opt] {base_name}: end", flush=True)

        # final diagnostics
        e_final = float(atoms.get_potential_energy())
        f_final = compute_fmax(atoms)
        converged = bool(f_final <= args.fmax)

        final_traj = outdir / f"{base_name}_final.traj"
        write(final_traj, atoms)

        result.update(
            {
                "status": "ok",
                "n_steps": int(getattr(opt, "nsteps", None) or 0),
                "e_final": e_final,
                "fmax_final": f_final,
                "converged": converged,
                "final_traj": str(final_traj),
                "logfile": str(logfile),
                "trajectory_log": str(trajlog),
            }
        )

    except Exception as exc:
        # Even on failure, try to capture whatever energy/forces we can
        result["status"] = "error"
        result["error"] = f"{type(exc).__name__}: {exc}"
        result["traceback"] = traceback.format_exc()

        try:
            # In case atoms exists and calc works, store partials
            if "atoms" in locals():
                result["e_partial"] = float(atoms.get_potential_energy())
                result["fmax_partial"] = compute_fmax(atoms)
        except Exception:
            pass

    finally:
        results_path.write_text(json.dumps(result, indent=2))
        print(f"[ok] wrote {results_path}", flush=True)


if __name__ == "__main__":
    main()
