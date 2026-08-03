#!/usr/bin/env python3
"""Triage xTB (GFN2) geometry optimisation runs for the Iridium-Oxide-ASE pipeline.

Walks a directory of run folders, parses the real artefacts produced by the
optimisation stage (ASE BFGS text logs, *_results.json, per site metadata JSON,
_opt.traj and _final.traj trajectories), and emits one metadata JSON per job
plus a single summary CSV across all jobs.

Outputs are written to a user specified out of repo path. This script is
committed to the repository; its outputs never are. Anyone cloning the repo
regenerates their own dataset and characterises it with this same instrument.

Triage bins:
    OK                  converged, no flags
    RESTART_MORE_STEPS  monotonic descent that hit the step ceiling,
                        restart from the LAST trajectory frame
    RESTART_FIRE        oscillation or plateau without progress,
                        switch BFGS to FIRE (soft adsorbate mode)
    RELABEL_SITE        H hopped sites, possibly a legitimate minimum,
                        the ML label must follow site_final
    DISCARD_PHYSICAL    desorption or surface reconstruction
    SETUP_ERROR         worst force on a deep slab atom, likely a missing
                        or inconsistent fixed layer constraint
    PENDING             structures and metadata exist but no optimiser
                        artefacts (job never ran)

SCF_FAIL is an annotation, not a bin: xTB electronic convergence warnings are
grepped from any xTB stdout files found beside the run. It rides alongside
whichever bin applies because its fix is electronic (temperature, damping),
not a different optimiser.

Usage:
    python triage_runs.py RUNS_DIR --out /path/outside/repo/triage_out \
        --config inputs/frozen_optimiser_config.json

Requires: ase, numpy, pandas.
"""

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import ase
    from ase.io import read as ase_read
    from ase.geometry import find_mic
    from ase.constraints import FixAtoms
except ImportError:
    sys.exit("ase is required: pip install ase")

SCRIPT_VERSION = "1.0.0"

# Window of optimiser steps inspected for descent, plateau and spike signals.
TAIL_WINDOW = 20
# Nearest surface neighbour beyond this distance, and increasing through the
# trajectory, reads as desorption.
DESORB_DISTANCE_A = 2.5
# Displacement of any free surface atom beyond this reads as reconstruction.
RECONSTRUCT_DISP_A = 1.5
# fmax jumping by more than this factor between adjacent steps is a spike.
SPIKE_FMAX_RATIO = 3.0
# Energy rising by more than this within the tail window is a spike (eV).
SPIKE_ENERGY_EV = 1.0
# Net fmax reduction over the tail window below this fraction is a plateau.
PLATEAU_NET_FRACTION = 0.10
# Fraction of decreasing steps in the tail window that counts as monotonic.
MONOTONIC_FRACTION = 0.60

LOG_ROW = re.compile(
    r"^\s*(?P<opt>[A-Za-z0-9_]+):?\s+(?P<step>\d+)\s+(?P<time>\d{2}:\d{2}:\d{2})"
    r"\s+(?P<energy>-?\d+\.\d+)\s+(?P<fmax>\d+\.\d+|\d+\.?\d*[eE][+-]?\d+)\s*$"
)

SCF_PATTERNS = re.compile(
    r"(scc.{0,40}not converged|scf.{0,40}not converged|convergence criteria"
    r".{0,40}not satisfied|no convergence|electronic temperature|warning)",
    re.IGNORECASE,
)


def sha256_short(path: Path, n: int = 12) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:n]


def parse_opt_log(log_path: Path) -> dict:
    """Parse an ASE optimiser text log into arrays plus wall time."""
    steps, energies, fmaxes, times, opt_names = [], [], [], [], []
    for line in log_path.read_text(errors="replace").splitlines():
        m = LOG_ROW.match(line)
        if m:
            opt_names.append(m.group("opt"))
            steps.append(int(m.group("step")))
            times.append(m.group("time"))
            energies.append(float(m.group("energy")))
            fmaxes.append(float(m.group("fmax")))
    wall_s = None
    if len(times) >= 2:
        # Timestamps are HH:MM:SS only; sum per contiguous segment and handle
        # midnight wrap. A restart appearing in the same log resets step
        # numbering, so segment on step decreases as well.
        total = timedelta()
        prev_t = datetime.strptime(times[0], "%H:%M:%S")
        for i in range(1, len(times)):
            t = datetime.strptime(times[i], "%H:%M:%S")
            if steps[i] <= steps[i - 1]:
                prev_t = t
                continue
            dt = t - prev_t
            if dt.total_seconds() < 0:
                dt += timedelta(days=1)
            total += dt
            prev_t = t
        wall_s = round(total.total_seconds(), 1)
    return {
        "optimiser": opt_names[-1] if opt_names else None,
        "steps": np.asarray(steps),
        "energies": np.asarray(energies),
        "fmaxes": np.asarray(fmaxes),
        "wall_time_s": wall_s,
    }


def fixed_indices(atoms) -> list:
    idx = []
    for c in atoms.constraints:
        if isinstance(c, FixAtoms):
            idx.extend(int(i) for i in c.index)
    return sorted(set(idx))


def mic_displacements(atoms_a, atoms_b) -> np.ndarray:
    """Per atom displacement between two frames, minimum image convention."""
    d = atoms_b.get_positions() - atoms_a.get_positions()
    vecs, _ = find_mic(d, atoms_a.get_cell(), pbc=atoms_a.get_pbc())
    return np.linalg.norm(vecs, axis=1)


def nearest_neighbour(atoms, index: int, exclude=()) -> dict:
    """Nearest atom to `index`, excluding given indices, MIC distances."""
    n = len(atoms)
    others = [i for i in range(n) if i != index and i not in exclude]
    d = atoms.get_distances(index, others, mic=True)
    j = others[int(np.argmin(d))]
    return {
        "index": int(j),
        "element": atoms[j].symbol,
        "distance_A": round(float(np.min(d)), 3),
    }


def grep_scf_warnings(job_dir: Path, base_name: str) -> dict:
    """Look for xTB stdout beside the run and grep electronic warnings.

    The current artefact bundles contain no xTB stdout, so absence is
    reported honestly rather than treated as a clean bill of health.
    """
    candidates = []
    for pattern in ("*.out", "*xtb*.log", "*stdout*", "*.txt"):
        candidates.extend(job_dir.glob(pattern))
    candidates = [c for c in candidates if base_name in c.name or "xtb" in c.name.lower()]
    if not candidates:
        return {"scf_output_found": False, "scf_warnings": []}
    warnings = []
    for c in candidates:
        for line in c.read_text(errors="replace").splitlines():
            if SCF_PATTERNS.search(line):
                warnings.append(f"{c.name}: {line.strip()[:160]}")
    return {"scf_output_found": True, "scf_warnings": warnings}


def tail_signals(energies: np.ndarray, fmaxes: np.ndarray) -> dict:
    """Signal 1: behaviour of fmax and energy over the last TAIL_WINDOW steps."""
    w = min(TAIL_WINDOW, len(fmaxes))
    f = fmaxes[-w:]
    e = energies[-w:]
    diffs = np.diff(f)
    frac_down = float(np.mean(diffs < 0)) if len(diffs) else 0.0
    net_fraction = float((f[0] - f[-1]) / f[0]) if f[0] > 0 else 0.0
    ratios = f[1:] / np.maximum(f[:-1], 1e-12)
    spike_fmax = bool(np.any(ratios > SPIKE_FMAX_RATIO))
    spike_energy = bool(np.any(np.diff(e) > SPIKE_ENERGY_EV))
    monotonic_descent = frac_down >= MONOTONIC_FRACTION and net_fraction > PLATEAU_NET_FRACTION
    plateau = net_fraction <= PLATEAU_NET_FRACTION and not (spike_fmax or spike_energy)
    return {
        "tail_window": int(w),
        "tail_frac_steps_decreasing": round(frac_down, 3),
        "tail_net_fmax_reduction_fraction": round(net_fraction, 3),
        "tail_spike_fmax": spike_fmax,
        "tail_spike_energy": spike_energy,
        "tail_monotonic_descent": bool(monotonic_descent),
        "tail_plateau_or_oscillation": bool(plateau),
    }


def analyse_job(results_path: Path, config_path: Path, config_hash: str) -> dict:
    job_dir = results_path.parent
    results = json.loads(results_path.read_text())
    base = results["base_name"]
    rec = {
        "job_id": base,
        "job_dir": str(job_dir),
        "config_hash": config_hash,
        "config_file": str(config_path),
        "ase_version": ase.__version__,
        "xtb_version": None,
        "script_version": SCRIPT_VERSION,
        "method": results.get("method"),
        "optimiser": None,
        "n_steps": results.get("n_steps"),
        "max_steps": results.get("max_steps"),
        "wall_time_s": None,
        "fmax_target": results.get("fmax_target"),
        "fmax_final": results.get("fmax_final"),
        "e_init": results.get("e_init"),
        "e_final": results.get("e_final"),
        "converged": bool(results.get("converged")),
        "status": "RUN",
    }
    cfg = json.loads(config_path.read_text())
    rec["xtb_version"] = cfg.get("xtb_version")

    # Own metadata only: each results dir carries metadata for every site,
    # so match on the job's base_name and never touch the others.
    meta_path = job_dir / f"metadata_{base}.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    # site_intended is preserved verbatim from the job's own metadata: the
    # raw value is carried unmodified and is never overwritten anywhere in
    # this script. The oNN label is derived alongside it for readability.
    rec["site_intended_raw"] = meta.get("target_o_index")
    rec["site_intended"] = (
        f"o{meta['target_o_index']}" if meta.get("target_o_index") is not None else None
    )
    rec["metadata_warnings"] = meta.get("warnings", [])
    z_freeze = meta.get("z_freeze_A")
    h_indices = meta.get("h_indices", [])

    # Optimiser log.
    log_path = job_dir / f"{base}_opt.log"
    log = parse_opt_log(log_path) if log_path.exists() else None
    if log:
        rec["optimiser"] = log["optimiser"]
        rec["wall_time_s"] = log["wall_time_s"]

    # Trajectory frames.
    traj_path = job_dir / f"{base}_opt.traj"
    final_path = job_dir / f"{base}_final.traj"
    first = last = None
    traj_len = 0
    if traj_path.exists():
        traj = ase_read(str(traj_path), index=":")
        traj_len = len(traj)
        first, last = traj[0], traj[-1]
    if final_path.exists():
        last = ase_read(str(final_path))

    # Constraint mask summary. The freeze rule (fix all atoms with z below
    # z_freeze_A) applies to the INITIAL structure: free atoms may relax
    # through the boundary later and that is physics, not a setup fault.
    if last is not None:
        fixed = fixed_indices(last)
        ref = first if first is not None else last
        z0 = ref.get_positions()[:, 2]
        rec["n_atoms"] = len(last)
        rec["n_fixed"] = len(fixed)
        rec["fixed_z_max_A"] = round(float(z0[fixed].max()), 3) if fixed else None
        rec["z_freeze_expected_A"] = z_freeze
        unfrozen_deep = (
            [int(i) for i in np.where(z0 < z_freeze)[0] if i not in fixed]
            if (z_freeze is not None)
            else []
        )
        rec["deep_atoms_not_fixed"] = unfrozen_deep[:20]
        rec["constraint_consistent_with_metadata"] = not unfrozen_deep
    else:
        fixed = []
        z0 = None

    # Worst force atom at the final frame, among free atoms, matching how the
    # optimiser computes fmax on constrained forces.
    rec["worst_force_atom"] = None
    forces_frame = None
    if last is not None:
        try:
            forces_frame = last.get_forces()
        except Exception:
            if traj_len:
                try:
                    forces_frame = traj[-1].get_forces()
                    last_for_forces = traj[-1]
                except Exception:
                    forces_frame = None
    if forces_frame is not None:
        fnorm = np.linalg.norm(forces_frame, axis=1)
        free = [i for i in range(len(fnorm)) if i not in fixed]
        wi = free[int(np.argmax(fnorm[free]))] if free else int(np.argmax(fnorm))
        rec["worst_force_atom"] = {
            "index": int(wi),
            "element": last[wi].symbol,
            "force_eV_per_A": round(float(fnorm[wi]), 4),
            "z_A": round(float(last.get_positions()[wi, 2]), 3),
            "z_initial_A": (round(float(z0[wi]), 3) if z0 is not None else None),
            "should_have_been_frozen": (
                bool(z0[wi] < z_freeze) if (z0 is not None and z_freeze is not None) else False
            ),
        }

    # Adsorbate tracking.
    rec["h_index"] = h_indices[0] if h_indices else None
    rec["h_displacement_A"] = None
    rec["max_free_surface_displacement_A"] = None
    rec["h_nn_start"] = rec["h_nn_end"] = None
    rec["site_final"] = None
    desorbing = False
    if first is not None and last is not None and len(first) == len(last):
        disp = mic_displacements(first, last)
        free_mask = np.ones(len(disp), bool)
        free_mask[fixed] = False
        for h in h_indices:
            free_mask[h] = False
        if free_mask.any():
            rec["max_free_surface_displacement_A"] = round(float(disp[free_mask].max()), 3)
        if h_indices:
            h = h_indices[0]
            rec["h_displacement_A"] = round(float(disp[h]), 3)
            rec["h_nn_start"] = nearest_neighbour(first, h, exclude=h_indices)
            rec["h_nn_end"] = nearest_neighbour(last, h, exclude=h_indices)
            rec["site_final"] = (
                f"o{rec['h_nn_end']['index']}"
                if rec["h_nn_end"]["element"] == "O"
                else f"{rec['h_nn_end']['element'].lower()}{rec['h_nn_end']['index']}"
            )
            # Desorption check: distance increasing through the trajectory.
            if traj_len >= 3:
                samples = np.linspace(0, traj_len - 1, min(5, traj_len)).astype(int)
                dists = []
                for s in samples:
                    fr = traj[s]
                    nn = nearest_neighbour(fr, h, exclude=h_indices)
                    dists.append(nn["distance_A"])
                rec["h_nn_distance_trace_A"] = dists
                desorbing = (
                    dists[-1] > DESORB_DISTANCE_A
                    and all(b >= a - 0.05 for a, b in zip(dists, dists[1:]))
                )

    # SCF warnings from any xTB stdout beside the run.
    rec.update(grep_scf_warnings(job_dir, base))
    scf_fail = bool(rec["scf_warnings"]) or any(
        "scf" in str(w).lower() or "scc" in str(w).lower() for w in rec["metadata_warnings"]
    )
    rec["scf_fail"] = scf_fail

    # Triage.
    signals = tail_signals(log["energies"], log["fmaxes"]) if log and len(log["fmaxes"]) else {}
    rec.update(signals)
    rec["triage_bin"], rec["triage_reason"] = triage(rec, signals, desorbing)
    if scf_fail:
        rec["triage_reason"] += "; SCF_FAIL annotation: fix is electronic (etemp or damping), not optimiser choice"
    return rec


def triage(rec: dict, signals: dict, desorbing: bool):
    hit_ceiling = (
        rec.get("n_steps") is not None
        and rec.get("max_steps") is not None
        and rec["n_steps"] >= rec["max_steps"]
    )
    wf = rec.get("worst_force_atom") or {}
    # Signal 2: an atom that started below the freeze boundary but is free,
    # and carries the worst force, means the fixed layer constraint is
    # missing or inconsistent. No optimiser choice fixes that.
    deep_worst = bool(wf.get("should_have_been_frozen"))
    missing_constraint = not rec.get("constraint_consistent_with_metadata", True)
    if deep_worst or missing_constraint:
        why = []
        if missing_constraint:
            why.append(
                f"atoms starting below z_freeze not fixed: {rec.get('deep_atoms_not_fixed')}"
            )
        if deep_worst:
            why.append(
                f"worst force on deep slab atom {wf.get('index')} ({wf.get('element')}), "
                f"initial z={wf.get('z_initial_A')}"
            )
        return "SETUP_ERROR", "; ".join(why)

    # Desorption is checked unconditionally: an H leaving the surface is not
    # rescued by any restart strategy.
    if desorbing:
        return "DISCARD_PHYSICAL", (
            f"H desorbing: nearest surface atom {rec['h_nn_end']['distance_A']} A and increasing"
        )

    hopped = (
        rec.get("site_intended")
        and rec.get("site_final")
        and rec["site_intended"] != rec["site_final"]
    )
    if hopped:
        return "RELABEL_SITE", (
            f"H moved {rec['site_intended']} -> {rec['site_final']}; possibly a legitimate "
            "minimum; ML label follows site_final, site_intended retained as metadata"
        )

    if rec.get("converged"):
        return "OK", "converged at intended site"

    max_disp = rec.get("max_free_surface_displacement_A")
    soft_mode_note = ""
    h = rec.get("h_index")
    if wf and h is not None:
        nn_end = rec.get("h_nn_end") or {}
        if wf.get("index") in (h, nn_end.get("index")):
            soft_mode_note = (
                f"; worst force on {'H' if wf.get('index') == h else 'H surface neighbour'} "
                "supports the soft mode reading"
            )

    # Signal 1 on the last steps of the optimiser log.
    if signals.get("tail_monotonic_descent") and hit_ceiling:
        return "RESTART_MORE_STEPS", (
            "monotonic descent hit the step ceiling; restart from the LAST trajectory frame, "
            "never from scratch"
        )
    if signals.get("tail_spike_fmax") or signals.get("tail_spike_energy"):
        # Spikes: consult signal 3, the geometry diff.
        if max_disp is not None and max_disp > RECONSTRUCT_DISP_A:
            return "DISCARD_PHYSICAL", (
                f"spiky forces with surface reconstruction: free surface atom moved {max_disp} A"
            )
        return "RESTART_FIRE", (
            "spikes in fmax or energy, geometry intact" + soft_mode_note + "; switch BFGS to FIRE"
        )
    if signals.get("tail_plateau_or_oscillation") or not signals.get("tail_monotonic_descent", False):
        note = ""
        if max_disp is not None and max_disp > RECONSTRUCT_DISP_A:
            note = f"; note large free surface atom displacement ({max_disp} A), review before restarting"
        return "RESTART_FIRE", (
            "oscillation or plateau without progress" + soft_mode_note + note + "; switch BFGS to FIRE"
        )
    return "RESTART_MORE_STEPS", "descending but unconverged; restart from last frame"


def discover_jobs(root: Path):
    """Find run jobs (results JSON) and pending jobs (metadata without results)."""
    results = sorted(root.rglob("*_results.json"))
    covered_dirs = {}
    for r in results:
        covered_dirs.setdefault(r.parent, set()).add(json.loads(r.read_text())["base_name"])
    pending = []
    for m in sorted(root.rglob("metadata_*.json")):
        base = m.stem.replace("metadata_", "")
        have = covered_dirs.get(m.parent, set())
        # Metadata bundles list every site; a metadata file is only a pending
        # job if no results exist for it anywhere under the root and the file
        # sits in a directory that owns no run for it. Restrict to metadata
        # whose own structure was the deliberate target: hydrated and pair
        # bundles carry exactly one such file per directory.
        if base not in have and not any(base in s for s in covered_dirs.values()):
            pending.append((m, base))
    # Deduplicate pending by base name, keep first occurrence.
    seen, pend = set(), []
    for m, base in pending:
        if base not in seen:
            seen.add(base)
            pend.append((m, base))
    return results, pend


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("runs_dir", type=Path, help="directory tree containing run folders")
    p.add_argument("--out", type=Path, required=True,
                   help="output directory OUTSIDE the repository; created if missing")
    p.add_argument("--config", type=Path, required=True,
                   help="frozen optimiser config (JSON), hashed into every record")
    p.add_argument("--csv-name", default="triage_summary.csv")
    args = p.parse_args()

    if not args.config.exists():
        sys.exit(f"config not found: {args.config}")
    cfg = json.loads(args.config.read_text())
    if str(cfg.get("xtb_version", "")).startswith("SET_ME"):
        print("WARNING: xtb_version in the frozen config is not pinned yet; "
              "fill it in from 'xtb --version' before the pair batch.", file=sys.stderr)
    config_hash = sha256_short(args.config)

    out = args.out.resolve()
    runs = args.runs_dir.resolve()
    if str(out).startswith(str(runs)):
        sys.exit("refusing to write outputs inside the runs directory")
    (out / "jobs").mkdir(parents=True, exist_ok=True)

    results, pending = discover_jobs(runs)
    if not results and not pending:
        sys.exit(f"no *_results.json or metadata_*.json found under {runs}")

    records = []
    for r in results:
        try:
            rec = analyse_job(r, args.config, config_hash)
        except Exception as e:
            rec = {"job_id": r.stem.replace("_results", ""), "job_dir": str(r.parent),
                   "status": "PARSE_ERROR", "triage_bin": "PARSE_ERROR",
                   "triage_reason": f"{type(e).__name__}: {e}", "config_hash": config_hash}
        records.append(rec)
    for m, base in pending:
        meta = json.loads(m.read_text())
        records.append({
            "job_id": base, "job_dir": str(m.parent), "status": "PENDING",
            "converged": None, "triage_bin": "PENDING",
            "triage_reason": "structures and metadata exist, never optimised",
            "site_intended": (f"o{meta['target_o_index']}"
                              if meta.get("target_o_index") is not None else None),
            "config_hash": config_hash,
        })

    # Per job JSON. Existing run directories are never touched.
    for rec in records:
        (out / "jobs" / f"{rec['job_id']}.json").write_text(json.dumps(rec, indent=2))

    # Constraint mask consistency across ALL jobs, converged included. The
    # comparison is on the fixed set itself (count and top z), not the total
    # atom count: a clean slab lacking the adsorbate H is expected to differ
    # in size while sharing the same frozen layers.
    masks = {}
    for rec in records:
        key = (rec.get("n_fixed"), rec.get("fixed_z_max_A"))
        if rec.get("n_fixed") is not None:
            masks.setdefault(key, []).append(rec["job_id"])
    if len(masks) > 1:
        print("\nCONSTRAINT MASK INCONSISTENCY across jobs:")
        for k, jobs in masks.items():
            print(f"  n_fixed={k[0]} fixed_z_max={k[1]}: {', '.join(jobs)}")
    elif masks:
        (k,) = masks.keys()
        print(f"\nConstraint masks consistent across all parsed jobs "
              f"(n_fixed={k[0]}, fixed_z_max={k[1]} A).")

    cols = ["job_id", "status", "site_intended", "site_final", "converged", "n_steps",
            "max_steps", "fmax_final", "fmax_target", "wall_time_s", "n_fixed",
            "constraint_consistent_with_metadata", "h_displacement_A", "scf_fail",
            "triage_bin", "triage_reason", "config_hash"]
    df = pd.DataFrame(records)
    for c in cols:
        if c not in df.columns:
            df[c] = None
    df = df[cols]
    csv_path = out / args.csv_name
    df.to_csv(csv_path, index=False)

    print(f"\n{len(records)} jobs -> {csv_path}")
    unconverged = df[(df["converged"] == False) & (df["status"] == "RUN")]  # noqa: E712
    if len(unconverged):
        print("\nTriage table, non converged jobs:")
        print(unconverged[["job_id", "site_intended", "site_final", "n_steps", "fmax_final",
                           "triage_bin", "triage_reason"]].to_string(index=False))
    print("\nReminder: restarts get NEW run ids; never modify an existing run directory.")


if __name__ == "__main__":
    main()
