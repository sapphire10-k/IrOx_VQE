#!/usr/bin/env python3
"""Render the boundary condition probe as a GitHub step summary, and fail the
build if the resolved environment contradicts the frozen optimiser config.

What is and is not a build failure:

  NOT a failure  GFN2-xTB refusing to run under periodic boundary conditions.
                 That is the result the probe exists to record. It exits 0 and
                 is reported as the verdict.
  NOT a failure  tblite segfaulting. Also a recorded result; the probe isolates
                 each trial in a subprocess precisely so this survives.
  FAILURE        the probe output is missing or malformed, which means the
                 probe did not run rather than that it found something.
  FAILURE        the installed xtb differs from xtb_version in the frozen
                 optimiser config. The config claims a frozen method; if the
                 environment has drifted from it the claim is false.

Usage:
    python scripts/ci_report.py \
        --probe outputs/pbc_feasibility_probe_linux.json \
        --conda-list conda-list.json \
        --config inputs/frozen_optimiser_config.json >> "$GITHUB_STEP_SUMMARY"
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PBC_LABEL = {"[true, true, true]": "TTT", "[true, true, false]": "TTF",
             "[false, false, false]": "FFF"}

# Packages worth naming in the summary. Anything else in the solve is noise.
OF_INTEREST = ("python", "numpy", "ase", "xtb", "xtb-python",
               "tblite", "tblite-python")


def pbc_label(pbc) -> str:
    return "".join("T" if p else "F" for p in pbc)


def load_conda_list(path: Path | None) -> dict:
    """conda list --json gives a list of package records."""
    if path is None or not path.exists():
        return {}
    data = json.loads(path.read_text())
    if isinstance(data, dict):           # some builds wrap it
        data = data.get("packages", data.get("result", []))
    return {rec["name"]: rec.get("version", "?")
            for rec in data if isinstance(rec, dict) and "name" in rec}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--probe", required=True)
    ap.add_argument("--conda-list", default=None)
    ap.add_argument("--config", default="inputs/frozen_optimiser_config.json")
    ap.add_argument("--vacuum", nargs="*", default=[],
                    help="Optional pbc_2d_vs_3d_check JSON files to summarise")
    args = ap.parse_args()

    out, failures = [], []

    probe_path = Path(args.probe)
    if not probe_path.exists():
        print(f"## Boundary condition probe\n\n"
              f"**Build failed.** The probe wrote no output at `{probe_path}`, "
              f"so it did not run to completion.")
        return 1
    try:
        probe = json.loads(probe_path.read_text())
        trials = probe["trials"]
    except (ValueError, KeyError) as exc:
        print(f"## Boundary condition probe\n\n"
              f"**Build failed.** `{probe_path}` is malformed: {exc}")
        return 1

    installed = load_conda_list(Path(args.conda_list) if args.conda_list else None)

    # The verdict: can GFN2-xTB be run periodically on this platform?
    deciding = [t for t in trials
                if t["calculator"] == "tblite" and t["method"] == "GFN2-xTB"
                and pbc_label(t["pbc"]) == "TTT"]
    out.append("## Can GFN2-xTB be run with periodic boundary conditions?\n")
    if deciding and deciding[0]["status"] == "ok":
        out.append("**Yes, through tblite.** The deciding trial "
                   "`tblite GFN2-xTB pbc=TTT` returned "
                   f"{deciding[0]['energy_eV']:.6f} eV.\n")
        out.append("Commits 10 and 11 in `docs/branch_commit_runbook.md` can "
                   "now be made. The Hamiltonian does not change, but three "
                   "settings do: the calculator moves from "
                   "`xtb.ase.calculator.XTB` to `tblite.ase.TBLite`, the "
                   "electronic temperature becomes an explicit recorded "
                   "argument instead of a library default, and "
                   "`vacuum_along_c_A` is set from the convergence tables "
                   "below rather than assumed. The pbc flag is not among "
                   "them: it was already set in commit 8.\n")
    elif deciding:
        d = deciding[0]
        detail = ("\n".join(d.get("error_chain") or [])
                  or d.get("error")
                  or f"{d.get('signal')} (returncode {d.get('returncode')})")
        out.append("**Not on this platform.** The deciding trial "
                   f"`tblite GFN2-xTB pbc=TTT` came back `{d['status']}`:\n")
        out.append(f"```\n{detail.strip()}\n```\n")
        out.append("Do not make commits 10 or 11. Commits 1 to 9 are correct "
                   "either way, so the branch still stands. See "
                   "`docs/pbc_recommendation.md` for the periodic DFT route.\n")
    else:
        out.append("The deciding trial was not present in the probe output.\n")

    # xtb refusal, quoted verbatim, because it is the whole finding.
    xtb_per = [t for t in trials
               if t["calculator"] == "xtb-python" and t["method"] == "GFN2-xTB"
               and pbc_label(t["pbc"]) == "TTT"]
    if xtb_per and xtb_per[0].get("status") != "ok":
        # Prefer the chained messages: the ASE wrapper's own string is a
        # generic "could not evaluate input" and says nothing useful.
        detail = "\n".join(xtb_per[0].get("error_chain")
                           or [xtb_per[0].get("error", "")])
        out.append("\n<details><summary>What xtb itself says about periodic "
                   "GFN2</summary>\n\n```\n" + detail.strip() + "\n```\n"
                   "</details>\n")

    # Environment, and the frozen config check.
    out.append("\n## Environment\n")
    out.append("| Package | Version |")
    out.append("| --- | --- |")
    for name in OF_INTEREST:
        if name in installed:
            out.append(f"| `{name}` | {installed[name]} |")
    if not installed:
        out.append("| (conda list not supplied) | |")

    cfg_path = Path(args.config)
    if cfg_path.exists() and "xtb" in installed:
        cfg = json.loads(cfg_path.read_text())
        want, got = str(cfg.get("xtb_version", "")), str(installed["xtb"])
        if want and want != got:
            failures.append(
                f"frozen config {cfg.get('config_name')} declares "
                f"xtb_version {want}, environment resolved xtb {got}")
            out.append(f"\n**Version drift.** `{cfg_path}` declares "
                       f"`xtb_version: {want}` but the environment resolved "
                       f"`xtb {got}`. The frozen method is not frozen.\n")
        elif want:
            out.append(f"\nFrozen config `{cfg.get('config_name')}` declares "
                       f"`xtb_version: {want}`, which matches the solve.\n")

    # Full trial matrix.
    out.append("\n## Trial matrix\n")
    out.append("| Engine | Method | pbc | Status | Energy (eV) | fmax (eV/A) |")
    out.append("| --- | --- | --- | --- | --- | --- |")
    for t in trials:
        e = t.get("energy_eV")
        f = t.get("fmax_all_atoms_eV_per_A")
        status = t["status"]
        if status == "crashed":
            status = f"crashed ({t.get('signal')})"
        out.append(
            f"| {t['calculator']} | {t['method']} | {pbc_label(t['pbc'])} | "
            f"{status} | {f'{e:.6f}' if e is not None else '-'} | "
            f"{f'{f:.4f}' if f is not None else '-'} |")

    # Optional vacuum convergence summaries.
    for vpath in args.vacuum:
        p = Path(vpath)
        if not p.exists():
            continue
        v = json.loads(p.read_text())
        out.append(f"\n## Vacuum sensitivity, `{p.name}`\n")
        out.append("| Engine | pbc | dE (eV) | Reading |")
        out.append("| --- | --- | --- | --- |")
        for s in v.get("summary", []):
            d = s.get("dE_for_vacuum_change_eV",
                      s.get("dE_for_20A_more_vacuum_eV"))
            out.append(f"| {s['calculator']} | {s['pbc']} | "
                       f"{f'{d:+.9f}' if d is not None else 'not determined'} | "
                       f"{s['verdict']} |")

    structure = probe.get("structure", {})
    if structure:
        out.append(f"\nSlab: {structure.get('formula')}, "
                   f"{structure.get('n_atoms')} atoms, vacuum along c "
                   f"{structure.get('vacuum_gap_along_c_A', float('nan')):.4f} A.\n")

    print("\n".join(out))

    if failures:
        print("\n".join(f"::error::{m}" for m in failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
