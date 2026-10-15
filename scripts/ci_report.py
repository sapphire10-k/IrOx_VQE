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
        --probe-conda-list conda-list-probe.json \
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
# dftd4 and multicharge are listed because GFN2 uses D4 dispersion and GFN1
# does not: a GFN2-only crash with GFN1 healthy points straight at them.
OF_INTEREST = ("python", "numpy", "ase", "xtb", "xtb-python",
               "tblite", "tblite-python", "dftd4", "multicharge")


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
    ap.add_argument("--probe-conda-list", default=None,
                    help="conda list --json for the throwaway probe env. Shown "
                         "beside the production env; never version-checked.")
    ap.add_argument("--config", default="inputs/frozen_optimiser_config.json")
    ap.add_argument("--vacuum", nargs="*", default=[],
                    help="Optional pbc_2d_vs_3d_check JSON files to summarise")
    args = ap.parse_args()

    out, failures = [], []

    probe_path = Path(args.probe)
    if not probe_path.exists():
        print(f"## Boundary condition probe\n\n"
              f"**No probe output** at `{probe_path}`. This step runs even "
              f"after an earlier failure, so this message usually means an "
              f"earlier step failed and the probe never started. The first "
              f"red step in the job log is the real error.")
        return 1
    try:
        probe = json.loads(probe_path.read_text())
        trials = probe["trials"]
    except (ValueError, KeyError) as exc:
        print(f"## Boundary condition probe\n\n"
              f"**Build failed.** `{probe_path}` is malformed: {exc}")
        return 1

    installed = load_conda_list(Path(args.conda_list) if args.conda_list else None)
    probe_env = load_conda_list(Path(args.probe_conda_list)
                                if args.probe_conda_list else None)

    # The verdict: can GFN2-xTB be run periodically on this platform?
    deciding = [t for t in trials
                if t["calculator"] == "tblite" and t["method"] == "GFN2-xTB"
                and pbc_label(t["pbc"]) == "TTT"]
    out.append("## Can GFN2-xTB be run with periodic boundary conditions?\n")
    if deciding and deciding[0]["status"] == "ok":
        out.append("**Yes, through tblite.** The deciding trial "
                   "`tblite GFN2-xTB pbc=TTT` returned "
                   f"{deciding[0]['energy_eV']:.6f} eV.\n")
        out.append("The engine change in "
                   "`docs/periodic_boundaries_change_set.md` can now be made. "
                   "The Hamiltonian does not change, but three settings do: "
                   "the calculator moves from `xtb.ase.calculator.XTB` to "
                   "`tblite.ase.TBLite`, the electronic temperature becomes "
                   "an explicit recorded argument instead of a library "
                   "default, and `vacuum_along_c_A` is set from the "
                   "convergence tables below rather than assumed. The pbc "
                   "flag is not among them: it is already in the frozen "
                   "config.\n")
    elif deciding:
        d = deciding[0]
        detail = ("\n".join(d.get("error_chain") or [])
                  or d.get("error")
                  or f"{d.get('signal')} (returncode {d.get('returncode')})")
        # The control that makes the periodic row interpretable: the same
        # engine and method with NO periodicity. If that also fails, the
        # periodic failure says nothing about periodicity; the GFN2 path in
        # this build is broken outright.
        control = [t for t in trials
                   if t["calculator"] == "tblite" and t["method"] == "GFN2-xTB"
                   and pbc_label(t["pbc"]) == "FFF"]
        control_failed = bool(control) and control[0]["status"] != "ok"
        if control_failed:
            out.append("**Not determined.** `tblite GFN2-xTB` failed with "
                       "periodic boundaries, but it also failed **without** "
                       f"them (`pbc=FFF` came back `{control[0]['status']}`). "
                       "A method that cannot run on an open cluster tells you "
                       "nothing about periodicity: the GFN2 path in this "
                       "tblite build is broken outright, and the question "
                       "this probe exists to answer has not been asked yet.\n")
            gfn1 = [t for t in trials
                    if t["calculator"] == "tblite" and t["method"] == "GFN1-xTB"
                    and pbc_label(t["pbc"]) == "TTT"]
            if gfn1 and gfn1[0]["status"] == "ok":
                out.append("tblite's periodic machinery itself works: "
                           "`tblite GFN1-xTB pbc=TTT` returned "
                           f"{gfn1[0]['energy_eV']:.6f} eV. The fault is "
                           "specific to the GFN2 path. GFN2 uses D4 "
                           "dispersion and GFN1 uses D3, so check the "
                           "resolved `dftd4` version in the probe column of "
                           "the environment table.\n")
            out.append(f"```\n{detail.strip()}\n```\n")
            out.append("Do not apply the engine change. Fix the build, then "
                       "re-run this workflow.\n")
        else:
            out.append("**No.** `tblite GFN2-xTB` runs on an open cluster but "
                       "fails with periodic boundaries, so the failure is "
                       f"about periodicity. It came back `{d['status']}`:\n")
            out.append(f"```\n{detail.strip()}\n```\n")
            out.append("Do not apply the engine change. See "
                       "`docs/pbc_recommendation.md` for the periodic DFT "
                       "route.\n")
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
    # Two environments, both built from environment.yml. iro2 is production
    # and is the only one version-checked against the frozen config. probe is
    # iro2 plus one override (the tblite pin), lives only inside the CI run,
    # and may legitimately resolve different versions, which is why the
    # differing rows are marked rather than failed.
    out.append("\n## Environment\n")
    if probe_env:
        out.append("| Package | `iro2` (production) | `probe` (this run's "
                   "calculations) |")
        out.append("| --- | --- | --- |")
        for name in OF_INTEREST:
            a, b = installed.get(name, "-"), probe_env.get(name, "-")
            if a == "-" and b == "-":
                continue
            flag = " **differs**" if a != b else ""
            out.append(f"| `{name}` | {a} | {b}{flag} |")
        out.append("\nThe trial matrix below was computed in `probe`. Only "
                   "`iro2` is checked against the frozen config.\n")
    else:
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
