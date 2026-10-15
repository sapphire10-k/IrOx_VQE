# Periodic boundary conditions: the change set, and how to land it

Single source of truth. This replaces `docs/pbc_code_changes.md` and
`docs/branch_commit_runbook.md`, which numbered the same work two different
ways and disagreed about the config edits. Delete both if they are already on
the branch.

Findings and alternatives live in `docs/pbc_recommendation.md`; measurements
in `docs/pbc_results_table.md`. This file is only the sequence.

**Division of labour.** You make and commit every change. Nothing here has
been applied to the repository, and no step below edits or commits anything
on your behalf.

**Execution status.** Steps 1 to 4 add files. Of the five new scripts, three
have been run and two have not; each carries an EXECUTION STATUS header, and
the "Known gaps" section at the end records exactly how far each was
exercised. Steps 5 to 11 are diffs to existing files, none of them applied.

---

## 0. Before you start

Create the branch: on the repository home page, click the branch dropdown
reading `main`, type `periodic-boundaries`, and click
**Create branch: periodic-boundaries from main**.

Three things to get right on every commit that follows.

- **Check the branch selector before each commit.** It must read
  `periodic-boundaries`. Arriving at a file through search or an old URL can
  put you back on `main` with no warning.
- **In the commit dialog choose "Commit directly to the periodic-boundaries
  branch".** Not the `main` option, and not the second radio, which would
  branch again off this one.
- **Replace the prefilled subject.** GitHub fills in "Update <filename>". The
  subject goes in the first box, the body in **Extended description**.

To edit a file: open it and press `e`, or click the pencil. To add files:
navigate **into** the target directory first, then **Add file, Upload files**.
Every file is already named as it should land; do not rename anything during
upload.

**Shortcut.** Press `.` with the repository open for github.dev, the browser
build of VS Code. It stages several files across several directories into one
commit, which collapses steps 1 to 4 into a single commit and makes the
directory boundaries irrelevant. The four-way split below exists only because
the plain upload form puts files in whichever single directory you are
standing in.

---

## 1. The sequence

| Step | File | Flow | Where | Conditional |
| --- | --- | --- | --- | --- |
| 1 | five new scripts | upload | `scripts/` | no |
| 2 | three docs, two result JSONs | upload | `docs/` | no |
| 3 | `environment.yml` | edit | repository root | no |
| 4 | `pbc-feasibility.yml` | upload | `.github/workflows/` | no |
| 5 | `generate_states.py` | edit | `scripts/` | no |
| 6 | `iro2_slab_setup.py` | edit | `scripts/` | no |
| 7 | `triage_surface_oxygens.py` | edit | `scripts/` | no |
| 8 | `frozen_optimiser_config.json` | edit | `inputs/` | no |
| 9 | `optimization_iro2_test.py` | edit | `scripts/` | no |
| | **read the verdict** | | Actions tab | |
| 10 | `frozen_optimiser_config.json` | edit | `inputs/` | **yes** |
| 11 | `optimization_iro2_test.py` | edit | `scripts/` | **yes** |
| | regenerate states, rerun relaxations | | not a commit | |

Steps 1 to 9 are correct whichever way the feasibility question goes, so the
branch stands on its own even if periodic GFN2 turns out to be unavailable.

Three ordering constraints, all load bearing in a one-file-at-a-time flow:

1. **Step 4 after steps 1 and 3.** The workflow runs `--help` on all five
   scripts and solves `environment.yml`, so both must exist when it first
   fires. Putting it fourth makes its first run also its first useful run.
2. **Step 8 before step 9.** The guard added in 9 reads `cfg["pbc"]`, which 8
   adds. A config field nothing reads breaks nothing; a reader for a field
   that does not exist raises `KeyError`.
3. **Step 10 before step 11**, for the same reason: 11 reads `--etemp` and
   the engine choice that 10 records.

**From step 8 the pipeline is deliberately inoperable**, and stays so until
step 11. After step 8 the config declares periodic GFN2-xTB, which xtb 6.7.1
cannot evaluate at all; after step 9 the guard additionally refuses every
structure currently on disk, all of which were built with periodicity off.
That is the intended state, not a regression. See the CI appendix.

---

## 2. Step 1. Five scripts to `scripts/`

```
pbc_feasibility_probe.py
pbc_2d_vs_3d_check.py
singlepoint_edge_vs_interior.py
relax_periodic.py
ci_report.py
```

> Add boundary-condition probes and the periodic relaxation scripts

```
pbc_feasibility_probe.py runs the 96-atom canonical clean slab through
xtb-python and tblite under (F,F,F), (T,T,F) and (T,T,T), for both GFN2 and
GFN1, and records the outcome of every combination. Each trial runs in its
own subprocess so that a library-level abort is recorded as a result rather
than taking the whole run down: tblite's GFN2 path segfaults on macOS arm64,
which an in-process loop would have hidden.

pbc_2d_vs_3d_check.py varies only the c lattice vector, with atoms held fixed
in Cartesian space, to distinguish genuine 2D periodicity from a third pbc
flag that is silently ignored.

singlepoint_edge_vs_interior.py compares open against periodic single points
on the unrelaxed canonical slab, reporting forces on the named edge O atoms
against interior O atoms under the frozen FixAtoms constraint, with raw
forces taken before the constraint zeroes the frozen rows. It uses the same
subprocess isolation and rewrites its output after every row, because a
tblite segfault is an OS-level crash no try/except can trap and would
otherwise discard the xtb-python rows already measured.

relax_periodic.py relaxes a structure under a chosen boundary condition with
the frozen optimiser configuration otherwise unchanged, and reports
convergence, minimum O-O and Ir-O, maximum mobile displacement, any O-O
contact below 1.60 A, and for an adsorbate run the H to site O distance,
nearest-O identity across the trajectory and a desorption verdict. All
distances use the minimum image convention wherever the corresponding pbc
flag is set. It never writes site_intended.

ci_report.py renders the probe output as a GitHub step summary. It has been
run against the macOS probe output, but only on the path where the probe
parses cleanly; its version drift branch, the one that can fail a build, has
not fired.

The four scripts that drive a calculator record the chained exception rather
than str(exc): ASE wraps a calculator failure in its own CalculationFailed,
so the library's own diagnosis, which is the entire finding here, sits one
level down in __context__. ci_report.py reads the recorded chain back out.

singlepoint_edge_vs_interior.py and relax_periodic.py have never been
executed. Both say so in their header.
```

## 3. Step 2. Five files to `docs/`

```
periodic_boundaries_change_set.md      (this file)
pbc_recommendation.md
pbc_results_table.md
pbc_feasibility_probe_macos_arm64.json
pbc_2d_vs_3d_check_macos_arm64.json
```

The two JSON files are the output of the first two scripts in step 1, not
copies of them. They carry the platform in the name because macOS arm64 is
exactly where tblite is broken and these numbers must not be mistaken for
Linux ones. The scripts' own `--out` defaults stay unsuffixed and point at
the untracked `outputs/`, so a fresh run never overwrites the committed
record.

> Document the periodic boundary findings and record the macOS measurements

```
Measured with xtb 6.7.1, xtb-python 22.1, tblite 0.6.0, ASE 3.29.0,
NumPy 2.4.6, Python 3.11.17 on macOS arm64.

GFN2-xTB cannot be run with periodic boundary conditions in xtb. The library
refuses with "scf: Multipoles not available with PBC": GFN2 carries atomic
dipoles and quadrupoles and that expansion has no periodic implementation.
GFN1, monopole only, runs periodically in both engines.

tblite 0.6.0 GFN2 segfaults on this platform even for non-periodic water
while its GFN1 works normally, so the tblite periodic-GFN2 question is
untested here and must be settled on Linux.

pbc (T,T,F) is not honoured as 2D by either engine: both return identical
energies for (T,T,F) and (T,T,T). Growing c by 20 A moves the xtb energy by
-0.7179 eV, so xtb uses the c vector and (T,T,F) is 3D in disguise; tblite
moves by -5.4e-7 eV and ignores c entirely.

The engines agree to 8.3e-6 eV on non-periodic GFN1 but disagree by 10.09 eV
on periodic GFN1 for the same structure.

The recommendation also records why keeping non-periodic relaxations behind a
cluster cut is not defensible for this cell. With a = 12.7231 A and
b = 6.3088 A and the boundaries open, no atom can be more than
min(a,b)/2 = 3.1544 A from a free face, against the pipeline's own 2.30 A
Ir-O neighbour criterion, and three of the four bridging sites lie exactly on
a boundary plane: o21 at (0.0000, 0.0000), o45 at (0.0000, 3.1544), o69 at
(6.3616, 0.0000). The cut is also applied after relaxation, so it removes the
broken atoms without removing their influence on where the survivors came to
rest.
```

## 4. Step 3. `environment.yml`, edit

```yaml
# before
  - xtb
  - xtb-python

# after
  - xtb=6.7.1          # matches xtb_version in the frozen optimiser config
  - xtb-python
  - tblite-python      # the candidate route for periodic GFN2
  - pyflakes           # catches undefined names; py_compile only parses
```

This triggers `python-package-conda.yml`, whose `generate-states` job builds
the environment and runs the existing pipeline. That is the right place to
discover whether the pin and the new dependency solve cleanly: before any
pipeline change, with an existing job as the test, when the only thing to
revert is one file. Only that one workflow fires; see the CI appendix.

> Pin xtb and add tblite to the shared environment

```
Three changes to the one environment every workflow already uses.

xtb is pinned to 6.7.1. inputs/frozen_optimiser_config.json declares
xtb_version 6.7.1 and triage_runs.py hashes that file into every metadata
record, but the environment resolved xtb unpinned, so the claim that the
method was frozen was not enforced anywhere. The feasibility workflow now
asserts the two agree on every run, and because that check runs against this
environment rather than a private one, it guards the environment the science
actually runs in.

tblite-python is added because xtb 6.7.1 cannot evaluate GFN2-xTB under
periodic boundary conditions at all, refusing with "scf: Multipoles not
available with PBC". tblite is the candidate route. It is left unpinned on
purpose: a release that adds periodic multipole support should surface in CI
as a changed version beside a changed verdict, not be held back by a pin.
Every run records the resolved versions, so the comparison stays auditable
without being frozen.

pyflakes is added because py_compile only parses. A name used solely on an
error path compiles cleanly and then raises NameError in the field, which is
a defect that got through review in this very branch.
```

## 5. Step 4. `pbc-feasibility.yml` to `.github/workflows/`

The directory already exists, so navigate into it and upload as before.

**This commit starts the run that answers the question.** If you would rather
it did not fire yet, delete the `push:` block before uploading and restore it
when you are ready; the workflow also accepts a manual dispatch.

> Add CI for the periodic boundary feasibility probe

```
GitHub's runners are Linux x86-64, which is the platform this question has to
be settled on: tblite's GFN2 path segfaults on macOS arm64 even for a water
molecule, so the local measurements cannot answer it.

This workflow is therefore not only a guard on the scripts. It is the Linux
test the investigation needs, and it reports the verdict to the job summary
where it can be read in the browser.

It runs in the shared environment.yml rather than a private one, so that the
xtb-version-against-frozen-config assertion guards the environment the
science actually runs in. A separate environment would have validated the
wrong thing.

scripts/ci_report.py fails the build on two conditions only: the probe
produced no parseable output, or the installed xtb differs from xtb_version
in the frozen optimiser config. A GFN2 refusal under periodic boundary
conditions is the finding, not a failure, and exits 0. So is a segfault,
which the probe's subprocess isolation turns into a recorded result.

The workflow also runs pyflakes over the five new scripts, because
py_compile only parses and a name used solely on an error path raises
NameError the first time that path is taken.
```

## 6. Step 5. `scripts/generate_states.py`, line 117

The states inherit the boundary condition every later stage uses, so this one
line is the origin of the artefact.

```python
# before
    # Read QE slab. We keep PBC off (matches your xTB/ASE pipeline style).
    slab0 = read(slab_in, format="espresso-in")
    slab0.set_pbc((False, False, False))

# after
    # Read QE slab. Periodic in a and b: the slab builder and the surface
    # triage gate both assume it, and open faces are only 3.15 A from the
    # cell centre, so nothing in this cell is free of edge effects.
    slab0 = read(slab_in, format="espresso-in")
    slab0.set_pbc((True, True, True))
```

Use `(True, True, True)`, not `(True, True, False)`. Neither engine honours a
mixed flag: both return bit identical energies for the two, and the xtb
result moves with the c axis, so `(T,T,F)` is 3D periodicity in disguise.

> Relax slabs periodically in a and b, not as finite clusters

```
generate_states.py set pbc (False, False, False) on the slab handed to every
downstream stage, while build_canonical_110_slab.py built it with
(True, True, False) and triage_surface_oxygens.py assumed the same. Nothing
compared them, so every relaxation treated the 96-atom slab as a finite
cluster with four free faces.

Uses (True, True, True) rather than (True, True, False): neither xtb nor
tblite honours a mixed flag, and the xtb result moves with the c axis, so
(True, True, True) records what is actually computed.

Existing relaxed structures were produced under the old boundary condition
and must be rebuilt, not relabelled. This commit regenerates nothing.

Second of two call sites follows in iro2_slab_setup.py.
```

## 7. Step 6. `scripts/iro2_slab_setup.py`, lines 60 and 171

Two more copies of the same call, not in your audit. Both need changing, or
states built through this path silently keep the old behaviour.

```python
# before
    slab.set_pbc((False, False, False))     # line 60
    slab0.set_pbc((False, False, False))    # line 171

# after
    slab.set_pbc((True, True, True))
    slab0.set_pbc((True, True, True))
```

> Fix the two remaining non-periodic set_pbc calls in iro2_slab_setup.py

```
Lines 60 and 171 carried further copies of the call corrected in
generate_states.py. States built through this path would otherwise have kept
the finite-cluster behaviour after the main fix landed. Not covered by the
original audit.
```

## 8. Step 7. `scripts/triage_surface_oxygens.py`, lines 55 to 58

Once states carry the right boundary condition, this restore becomes a silent
no-op that hides a mismatch rather than reporting one.

```python
# before
    # slabs with pbc switched off (xTB convention), so restore x/y
    atoms.set_pbc((True, True, False))

# after
    # States are written periodic in a and b. Assert rather than restore, so
    # a structure built by an older path is rejected instead of patched.
    assert atoms.get_pbc()[0] and atoms.get_pbc()[1], (
        f"expected a,b periodic, got pbc={atoms.get_pbc()}")
```

> Assert a,b periodicity in the surface-oxygen gate instead of restoring it

```
The gate silently reset pbc to (True, True, False) to compensate for states
written with periodicity switched off. With states now built periodic that
line is a no-op which would hide a mismatch rather than report one, so it
becomes an assertion. A structure built by an older path is now rejected
instead of quietly patched.
```

## 9. Step 8. `inputs/frozen_optimiser_config.json`, edit

Two keys change and one is added. **Do not add the engine or the vacuum
here**; those are step 10 and depend on the verdict.

```json
// before
  "config_name": "irox_gfn2_bfgs_v3",
  "description": "Frozen optimiser settings for all IrO2(110) H adsorption runs. ...",

// after
  "config_name": "irox_bfgs_v4_periodic_pending_engine",
  "description": "Frozen optimiser settings, periodic in a and b. Supersedes irox_gfn2_bfgs_v3, which relaxed every slab as a finite cluster. The engine is not yet settled: xtb 6.7.1 cannot evaluate GFN2-xTB under PBC.",
  "pbc": [true, true, true],
```

Everything else stays as it is, `xtb_version` included. The name says
`pending_engine` on purpose: until the probe reports, this config declares a
method that no installed engine can execute, and the name should not pretend
otherwise. `triage_runs.py` hashes this file into every metadata record, so
the rename is what keeps the old runs distinguishable under
`irox_gfn2_bfgs_v3`.

> Record the boundary condition in the frozen optimiser config

```
The boundary condition is part of the frozen method, so it belongs in the
config rather than being implied by whichever script built the structure.
triage_runs.py hashes this file into every metadata record, so bumping
config_name is what keeps the old runs distinguishable under
irox_gfn2_bfgs_v3.

gfn_level and xtb_version are deliberately unchanged and no engine key is
added. Whether GFN2 can be run periodically at all depends on the tblite
result, which is not yet in, so config_name carries pending_engine.

The guard that enforces this field follows in optimization_iro2_test.py.
```

## 10. Step 9. `scripts/optimization_iro2_test.py`, add the guard

Add a `--config` argument:

```python
    ap.add_argument("--config", default="inputs/frozen_optimiser_config.json",
                    help="Frozen optimiser config, checked against the structure")
```

and, immediately after `atoms = read(args.traj)`:

```python
        cfg = json.loads(Path(args.config).read_text())
        expected = tuple(cfg["pbc"])
        got = tuple(bool(p) for p in atoms.get_pbc())
        if got != expected:
            raise SystemExit(
                f"{args.traj}: pbc {got} does not match frozen config "
                f"{cfg['config_name']} which requires {expected}. "
                f"Rebuild the states with the current generate_states.py.")
```

`json` and `Path` are already imported in that file, so no new imports.

> Refuse to optimise a structure whose pbc contradicts the frozen config

```
The original fault was three files disagreeing about the boundary condition
with nothing comparing them. This compares the structure's pbc against the
config at the point structures enter the optimiser and exits on a mismatch,
so that class of bug is loud rather than silent.

Adds a --config argument defaulting to inputs/frozen_optimiser_config.json.
Requires the pbc field added in the preceding commit.

From this commit onwards the pipeline refuses every structure currently on
disk, all of which were built with periodicity off. Regenerate the states
before running anything.
```

---

## 11. Read the verdict

Open the pull request now if you have not already: accept the
"Compare and pull request" banner and open it as a **draft**.

> Periodic boundary conditions for the slab relaxations

```
Every relaxation in the pipeline has been treating the 96-atom slab as a
finite cluster with four free faces, which is the source of the 1.31 to
1.38 A O-O contacts and the 1.57 A Ir-O contacts in the audit.

Draft until the feasibility workflow has reported. That run decides whether
periodic GFN2-xTB is available at all, and therefore whether the last two
commits on this branch can be made.

See docs/pbc_recommendation.md for the findings and the alternatives.
```

You do not need a Linux machine. Step 4 already started the run. Open the
**Actions** tab, select **Periodic boundary feasibility (GFN2-xTB)**, and
read the job summary: the verdict in plain words, then the resolved package
versions, the trial matrix and the vacuum sensitivity tables. The raw JSON is
attached as the `pbc-probe-linux` artifact for 90 days. If no run appears,
use **Run workflow**.

The deciding row is `tblite GFN2-xTB pbc=TTT`.

- **`ok`**: periodic GFN2 exists. Continue to steps 10 and 11.
- **anything else**: periodic GFN2 is unavailable in the xtb and tblite
  family. **Stop here.** Do not make steps 10 or 11. Steps 1 to 9 are correct
  regardless, so the branch still stands, and the decision moves to the
  periodic DFT route in `docs/pbc_recommendation.md`. Note that the pipeline
  remains inoperable in that state, which is the honest outcome rather than
  something to patch around.

---

## 12. Step 10. `inputs/frozen_optimiser_config.json` again, conditional

Take `vacuum_along_c_A` from the GFN2 vacuum table in the job summary rather
than assuming the canonical 30.000 A is converged; xtb's GFN1 energy moved
0.718 eV for 20 A more vacuum. If that table reads "not determined", the step
was killed by a segfault and you will need to widen the sweep by hand with
`--c0` and `--c1`.

```json
// changed
  "config_name": "irox_gfn2_bfgs_v4_periodic",
  "description": "Frozen optimiser settings, periodic, tblite engine. Supersedes irox_gfn2_bfgs_v3.",

// added
  "engine": "tblite",
  "tblite_version": "<as resolved in the job summary>",
  "vacuum_along_c_A": <from the GFN2 vacuum table>,
  "k_points": "gamma",
  "electronic_temperature_K": 300.0,
  "accuracy": 1.0,
```

`k_points` records gamma because that is what tblite 0.6.0 offers. If you
need sampling, that is an argument for DFTB+ or CP2K rather than this route.

> Set the periodic engine and converged vacuum in the frozen config

```
Conditional on the feasibility workflow reporting tblite GFN2-xTB pbc=TTT
as ok.

vacuum_along_c_A comes from the convergence test rather than assuming the
canonical 30.000 A is converged: xtb GFN1 moved 0.718 eV for 20 A more
vacuum. k_points records gamma because that is what tblite 0.6.0 offers.

config_name drops pending_engine now that the engine is settled.
```

## 13. Step 11. `scripts/optimization_iro2_test.py` again, conditional

```python
# line 11, before
from xtb.ase.calculator import XTB
# after
from tblite.ase import TBLite

# the calculator line, before
        atoms.calc = XTB(method=args.method)
# after
        atoms.calc = TBLite(
            method=args.method,
            electronic_temperature=args.etemp,
            accuracy=1.0,
        )
```

plus the argument, so the electronic temperature is recorded rather than
inherited:

```python
    ap.add_argument("--etemp", type=float, default=300.0,
                    help="Electronic temperature in K, recorded in results JSON")
```

> Switch the optimiser to tblite for periodic GFN2-xTB

```
Conditional on the feasibility workflow reporting tblite GFN2-xTB pbc=TTT
as ok.

xtb 6.7.1 cannot evaluate GFN2 under periodic boundary conditions at all, so
the calculator moves from xtb.ase.calculator.XTB to tblite.ase.TBLite.

Electronic temperature becomes an explicit argument recorded in the results
JSON rather than inherited from whichever library is in use. The two codes
agree to 8.3e-6 eV on non-periodic GFN1, which indicates their defaults
currently match, but that is not something to depend on silently.
```

## 14. After step 11: regenerate, do not patch

Not a commit. Existing relaxed structures were produced under the old
boundary condition and cannot be repaired by relabelling. Regenerate the
states and rerun the relaxations. Do not overwrite `site_intended` at any
point, and keep the old outputs under the `irox_gfn2_bfgs_v3` config name so
the two sets stay distinguishable. `outputs/` is untracked, so record the
regeneration in `REBUILD_TRACKER.md` rather than committing files.

Then steps 2, 3 and 4 of the original investigation become runnable:
`singlepoint_edge_vs_interior.py` for the edge-versus-interior forces, and
`relax_periodic.py` for the clean slab and, only if that comes back clean,
for H*-o69.

---

## Known gaps in what you are committing

- The failure path in `pbc_2d_vs_3d_check.py` is **not exercised locally**.
  Its success path reproduces the committed numbers exactly and the
  `--calculators` flag is used by the workflow, but the branch that records a
  refusal rather than raising is first exercised by the GFN2 vacuum step in
  CI. That step carries `|| true` for the same reason: the script has no
  subprocess isolation, so a tblite segfault would otherwise kill it.
- `singlepoint_edge_vs_interior.py` and `relax_periodic.py` have never been
  run. CI checks that they compile, that pyflakes finds no undefined names,
  and that their argument parsers work, nothing more. In particular the
  subprocess isolation in `singlepoint_edge_vs_interior.py` has not been
  exercised against an actual crash.
- `ci_report.py` has been exercised only against the macOS probe output. Its
  version drift branch, the one that can fail a build, has not fired.
- The edge O list in `singlepoint_edge_vs_interior.py` (65, 70, 82, 89, 14,
  16, 56, 63) is taken from the audit on the assumption that it uses the same
  zero based ASE indexing as the site labels, which was confirmed against
  `state_set.json`. The list itself was not independently confirmed. The
  script prints `d_edge_A` per named atom, so check those are small on first
  run; anything near 3.15 A means the convention differs.

## Appendix: what this does to the ten existing workflows

Less than the path filters suggest, because two of the relevant workflows
restrict pushes by branch. Computed against the actual `on:` blocks:

| Step | Fires on push to this branch |
| --- | --- |
| 1, 2 | nothing |
| 3 | `python-package-conda.yml` |
| 4 | `pbc-feasibility.yml` |
| 5 | `python-package-conda.yml` |
| 6, 7 | nothing |
| 8 | `pbc-feasibility.yml` |
| 9 | `python-package-conda.yml` |

Opening the pull request additionally fires `example_run_all.yml` and
`pbc-feasibility.yml`, both of which carry a `pull_request` trigger.

**`example_run_all.yml` will go red, and it should.** It restricts pushes to
`main`, so it never fires while you commit, but its `pull_request` trigger
watches `scripts/**` and it runs the whole canonical example end to end:
`iro2_slab_setup.py`, then `optimization_iro2_test.py`, then cluster
extraction and orbital analysis. From step 6 onwards `iro2_slab_setup.py`
writes periodic structures, xtb cannot evaluate GFN2 under periodic boundary
conditions, the optimiser writes no structure, and
`examples/iro2_h_star_o69/run_all.sh` exits 1 on its own
"No optimised .traj found" check. That script does use `|| true` in six live
places, but every one of them is on a listing or filename-discovery command
(`ls -la`, `find | head`, `ls | head`); none is on a calculation, and none
guards the `exit 1`. Expect this PR check to be red from the moment you open
it until step 11 lands, and read it as the pipeline telling the truth.

**`python-package-conda.yml` has the opposite failure mode.** It triggers on
pushes touching `scripts/generate_states.py`, `environment.yml` and
`scripts/optimization_iro2_test.py`, so steps 3, 5 and 9 start it. Its
`generate-states` job should pass: the gate runs `triage_surface_oxygens.py`
against the `.in` file, and ASE reads a Quantum ESPRESSO input as fully
periodic, so the assertion added in step 7 is satisfied either way.

The trap is in its other two jobs. `optimize-single` and `optimize-clean-ref`
run `optimization_iro2_test.py` with `xtb.ase.calculator.XTB` and GFN2. Both
are `workflow_dispatch` only, so nothing happens automatically, but both wrap
the optimisation in `|| true`, which means that if you dispatch them they go
**green while producing results JSON whose status is "error"**. Green CI, no
results. Do not dispatch those two on this branch until step 11 lands.

`pseudopotential-qe-inputs.yml` watches `environment.yml` and
`generate_states.py` but restricts pushes to
`dev/pseudopotential-integration`, so it never fires here. It runs no xTB
optimisation in any case.

Counting up: of the ten pre-existing workflows, one fires on push
(`python-package-conda`), one when the pull request opens
(`example_run_all`), and `pseudopotential-qe-inputs` is branch-restricted
away. The remaining seven fire for nothing this branch touches:
`cluster-extraction`, `exact_diagonalization`, `export_hamiltonian_manual`,
`hamiltonian-export`, `orbital-analysis` and `triage` are
`workflow_dispatch` only, and `draft-pdf` triggers on push but watches
`paper.md`, `paper.bib` and its own file, none of which change here.
