# Browser runbook: landing the periodic boundary work on `periodic-boundaries`

Everything below is done on github.com. Nothing here commits anything for you.
Work top to bottom; the ordering is load bearing in two places, flagged where
it matters.

## 0. Create the branch

On the repository home page, click the branch dropdown that currently reads
`main`, type `periodic-boundaries` into the box, and click
**Create branch: periodic-boundaries from main**. You are now on the branch.

## 1. Browser mechanics, once

Three things to get right on every commit that follows.

- **Check the branch selector before each commit.** It must read
  `periodic-boundaries`. Arriving at a file through search or an old URL can
  put you back on `main` without any warning.
- **In the commit dialog, choose "Commit directly to the periodic-boundaries
  branch".** Never the `main` option, and not the second radio either, which
  would start a further branch off this one.
- **Replace the prefilled subject.** GitHub fills in "Update <filename>". The
  subject goes in the first box, the body in **Extended description**.

To edit an existing file: open it and press `e`, or click the pencil.
To add files: navigate **into** the target directory first, then
**Add file, Upload files**.

Every file below is already named as it should land. Do not rename anything
during upload.

### The one shortcut worth knowing

Press `.` with the repository open and you get github.dev, the browser build
of VS Code. It stages several files across several directories into a single
commit. If you use it, commits 1 to 4 below collapse into one and the
directory boundaries stop mattering. The split into four exists only because
the plain upload form puts files in whichever single directory you are
standing in.

## 2. Commit order

| # | Files | Flow | Where |
| --- | --- | --- | --- |
| 1 | the five new scripts | upload | `scripts/` |
| 2 | four markdown docs, two result JSONs | upload | `docs/` |
| 3 | `environment-pbc.yml` | upload | repository root |
| 4 | `pbc-feasibility.yml` | upload | `.github/workflows/` |
| 5 | `generate_states.py` | edit | `scripts/` |
| 6 | `iro2_slab_setup.py` | edit | `scripts/` |
| 7 | `triage_surface_oxygens.py` | edit | `scripts/` |
| 8 | `frozen_optimiser_config.json` | edit | `inputs/` |
| 9 | `optimization_iro2_test.py` | edit | `scripts/` |

Commit 4 starts the run that answers the whole question. Then stop: commits
10 and 11 depend on what it says.

**Ordering constraint 1: 8 before 9.** The guard added in 9 reads
`cfg["pbc"]`, which 8 adds. A config field nothing reads breaks nothing; a
reader for a field that does not exist raises `KeyError`.

**Ordering constraint 2: from 9 onwards the branch will not run against your
existing states.** The guard rejects any structure whose pbc disagrees with
the config, and every structure currently on disk was built with periodicity
off. That is the intended behaviour, but do not be surprised by it.

**Why 1 to 4 are in that order.** The workflow runs `--help` on every script
and solves the environment file, so both must already exist when it first
fires. Putting the workflow last means its first run is also its first
meaningful run.

## 3. The commits

### 1. `scripts/`, upload five files

```
pbc_feasibility_probe.py
pbc_2d_vs_3d_check.py
singlepoint_edge_vs_interior.py
relax_periodic.py
ci_report.py
```

Three of these have been run and two have not. That distinction lives in each
file's EXECUTION STATUS header, which is where it belongs, rather than in the
commit history. Section 6 records exactly how far each was exercised.

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
subprocess isolation as the probe, and rewrites its output after every row,
because a tblite segfault is an OS-level crash that no try/except can trap
and would otherwise discard the xtb-python rows already measured.

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

### 2. `docs/`, upload six files

```
pbc_recommendation.md
pbc_results_table.md
pbc_code_changes.md
branch_commit_runbook.md
pbc_feasibility_probe_macos_arm64.json
pbc_2d_vs_3d_check_macos_arm64.json
```

The two JSON files are the output of the first two scripts in commit 1, not
copies of them. They carry the platform in the name because macOS arm64 is
exactly where tblite is broken, and these numbers must not be mistaken for
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

### 3. Repository root, upload one file

`environment-pbc.yml`.

> Add a pinned environment for the boundary condition investigation

```
Deliberately separate from environment.yml. The production environment is
shared by ten existing workflows, and adding tblite or pinning xtb there
would change the behaviour of all of them.

Versions are pinned to those used for the macOS measurements, so that a
difference between a CI run and docs/pbc_feasibility_probe_macos_arm64.json
is a platform difference rather than a version difference.
```

### 4. `.github/workflows/`, upload one file

`pbc-feasibility.yml`. The directory already exists, so navigate into it and
use the upload form as before.

**This commit starts the run that answers the question.** If you would rather
it did not fire yet, delete the `push:` block before uploading and add it
back when you are ready; the workflow also accepts a manual dispatch.

> Add CI for the periodic boundary feasibility probe

```
GitHub's runners are Linux x86-64, which is the platform this question has to
be settled on: tblite's GFN2 path segfaults on macOS arm64 even for a water
molecule, so the local measurements cannot answer it.

This workflow is therefore not only a guard on the scripts. It is the Linux
test that step 1 of docs/branch_commit_runbook.md calls for, and it reports
the verdict to the job summary where it can be read in the browser.

scripts/ci_report.py fails the build on two conditions only: the probe
produced no parseable output, or the installed xtb differs from xtb_version
in the frozen optimiser config. A GFN2 refusal under periodic boundary
conditions is the finding, not a failure, and exits 0. So is a segfault,
which the probe's subprocess isolation turns into a recorded result.

The workflow also runs pyflakes over the five new scripts. py_compile only
parses, so a name used solely on an error path, such as sys.executable in a
crash handler, compiles cleanly and then raises NameError the first time
that path is taken. That exact defect was shipped and caught during review.
```

### 5. `scripts/generate_states.py`, line 117

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

### 6. `scripts/iro2_slab_setup.py`, lines 60 and 171

> Fix the two remaining non-periodic set_pbc calls in iro2_slab_setup.py

```
Lines 60 and 171 carried further copies of the call corrected in
generate_states.py. States built through this path would otherwise have kept
the finite-cluster behaviour after the main fix landed. Not covered by the
original audit.
```

### 7. `scripts/triage_surface_oxygens.py`, lines 55 to 58

> Assert a,b periodicity in the surface-oxygen gate instead of restoring it

```
The gate silently reset pbc to (True, True, False) to compensate for states
written with periodicity switched off. With states now built periodic that
line is a no-op which would hide a mismatch rather than report one, so it
becomes an assertion. A structure built by an older path is now rejected
instead of quietly patched.
```

### 8. `inputs/frozen_optimiser_config.json`

> Record the boundary condition in the frozen optimiser config

```
The boundary condition is part of the frozen method, so it belongs in the
config rather than being implied by whichever script built the structure.
triage_runs.py hashes this file into every metadata record, so bumping
config_name is what keeps the old runs distinguishable under
irox_gfn2_bfgs_v3.

gfn_level and engine are deliberately left unchanged. Whether GFN2 can be run
periodically at all depends on the tblite result, which is not yet in.

The guard that enforces this field follows in optimization_iro2_test.py.
```

### 9. `scripts/optimization_iro2_test.py`

Adds the `--config` argument as well as the guard.

> Refuse to optimise a structure whose pbc contradicts the frozen config

```
The original fault was three files disagreeing about the boundary condition
with nothing comparing them. This compares the structure's pbc against the
config at the point structures enter the optimiser and exits on a mismatch,
so that class of bug is loud rather than silent.

Adds a --config argument defaulting to inputs/frozen_optimiser_config.json.
Requires the pbc field added in the preceding commit.

From this commit onwards the pipeline will refuse every structure currently
on disk, all of which were built with periodicity off. Regenerate the states
before running anything.
```

## 4. Stop here and read the answer

Go to the repository, accept the "Compare and pull request" banner, and open
it as a **draft**. Suggested title and body:

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

You do not need a Linux machine. Commit 4 already started the run.

Open the **Actions** tab, select **Periodic boundary feasibility (GFN2-xTB)**,
and read the job summary. It opens with the verdict in plain words, followed
by the resolved package versions, the full trial matrix and the vacuum
sensitivity tables. The raw JSON is attached to the run as the
`pbc-probe-linux` artifact, kept for 90 days.

If no run appears, use **Run workflow** on that page.

The deciding row is `tblite GFN2-xTB pbc=TTT`.

- **`ok`**: periodic GFN2 exists. Continue to commits 10 and 11.
- **anything else**: periodic GFN2 is unavailable in the xtb and tblite
  family. Do not make commits 10 or 11. The branch still stands on its own,
  since commits 1 to 9 are correct either way, and the conversation moves to
  the periodic DFT route in `docs/pbc_recommendation.md`.

## 5. Conditional commits

Only if the workflow said `ok`. Same ordering rule, config before code.

### 10. `inputs/frozen_optimiser_config.json`

The vacuum convergence has already been measured for you: the workflow runs
`pbc_2d_vs_3d_check.py` under both GFN1 and GFN2 and puts both tables in the
job summary. Take `vacuum_along_c_A` from the GFN2 table rather than assuming
the canonical 30.000 A is converged. If the GFN2 table reads "not determined",
that step was killed by a segfault and you will need to widen the c sweep
manually with `--c0` and `--c1`.

> Set the periodic engine and converged vacuum in the frozen config

```
Conditional on the feasibility workflow reporting tblite GFN2-xTB pbc=TTT
as ok.

vacuum_along_c_A comes from the convergence test rather than assuming the
canonical 30.000 A is converged: xtb GFN1 moved 0.718 eV for 20 A more
vacuum. k_points records gamma because that is what tblite 0.6.0 offers.
```

### 11. `scripts/optimization_iro2_test.py`

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

## 6. Known gaps in what you are committing

- The failure path in `pbc_2d_vs_3d_check.py` is **not exercised locally**.
  Its success path reproduces the committed numbers exactly, and the
  `--calculators` flag is used by the workflow, but the branch that records a
  refusal rather than raising is first exercised by the GFN2 vacuum step in
  CI. That step carries `|| true` for the same reason: the script has no
  subprocess isolation, so a tblite segfault would otherwise kill it.
- `singlepoint_edge_vs_interior.py` and `relax_periodic.py` have never been
  run at all. CI checks that they compile, that pyflakes finds no undefined
  names, and that their argument parsers work, nothing more. In particular
  the subprocess isolation in `singlepoint_edge_vs_interior.py` has not been
  exercised against an actual crash.
- `ci_report.py` has been exercised only against the macOS probe output. Its
  version drift branch, which is the one that can fail a build, has not
  fired.
- The edge O list in `singlepoint_edge_vs_interior.py` (65, 70, 82, 89, 14,
  16, 56, 63) is taken from the audit on the assumption that it uses the same
  zero based ASE indexing as the site labels, which was confirmed against
  `state_set.json`. The list itself was not independently confirmed. The
  script prints `d_edge_A` per named atom, so check those are small on first
  run; anything near 3.15 A means the convention differs.
- State regeneration is not a commit here, since `outputs/` is untracked.
  Record it in `REBUILD_TRACKER.md` instead.

## 7. What this branch does to the ten existing workflows

Only one of them reacts, and it needs watching.

`python-package-conda.yml` triggers on pushes touching
`scripts/generate_states.py`, so **commit 5 will start it**. Its
`generate-states` job should pass: the gate runs
`triage_surface_oxygens.py` against the `.in` file, and ASE reads a Quantum
ESPRESSO input as fully periodic, so the assertion added in commit 7 is
satisfied either way.

The trap is in its other two jobs. `optimize-single` and
`optimize-clean-ref` run `optimization_iro2_test.py`, which uses
`xtb.ase.calculator.XTB` with GFN2. From commit 5 onwards the states they
optimise are periodic, and **xtb cannot evaluate GFN2 under periodic boundary
conditions at all**. Every optimisation in those jobs will fail.

Both jobs only run on `workflow_dispatch`, so nothing happens automatically.
But both wrap the optimisation in `|| true`, which means that if you do
dispatch them the jobs will go **green while producing results JSON files
whose status is "error"**. Green CI, no results. Do not dispatch those two
jobs on this branch until the engine question is settled, and when it is,
`optimization_iro2_test.py` needs the commit 11 calculator swap before they
mean anything again.

The other nine workflows do not touch the changed files and are unaffected.
