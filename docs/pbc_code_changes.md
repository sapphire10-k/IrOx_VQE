# Exact code changes implied by the recommendation

You make and commit these yourself. Nothing below has been applied to the
repository. Changes 1 to 3 and 6 are worth making whatever the tblite answer
turns out to be; changes 4 and 5 depend on it.

## 0. Settle this first

```bash
python scripts/pbc_feasibility_probe.py          # on Linux, not macOS
```

If `tblite GFN2-xTB pbc=TTT` reports `ok`, periodic GFN2 is available and the
Hamiltonian does not change. If it reports a multipole error like xtb's, it
is not, and the recommendation moves to periodic DFT.

## 1. `scripts/generate_states.py`, line 117

The states inherit the boundary condition that every later stage uses, so
this single line is the origin of the artefact.

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

Use `(True, True, True)`, not `(True, True, False)`. Neither xtb nor tblite
honours a mixed flag: both return bit identical energies for the two, and the
xtb result moves with the c axis, so `(T,T,F)` is 3D periodicity in disguise.
Writing `(T,T,T)` records what is actually computed.

## 2. `scripts/iro2_slab_setup.py`, lines 60 and 171

Two more copies of the same call, not mentioned in your audit. Both need the
same change, or states built through this path will silently keep the old
behaviour.

```python
# both sites, before
    slab.set_pbc((False, False, False))     # line 60
    slab0.set_pbc((False, False, False))    # line 171

# after
    slab.set_pbc((True, True, True))
    slab0.set_pbc((True, True, True))
```

## 3. `scripts/triage_surface_oxygens.py`, lines 55 to 58

Once states carry the correct boundary condition this restore becomes a
silent no-op that hides a mismatch rather than reporting one.

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

## 4. `scripts/optimization_iro2_test.py`, lines 11 and 52

Only if the probe shows tblite gives periodic GFN2.

```python
# line 11, before
from xtb.ase.calculator import XTB
# after
from tblite.ase import TBLite

# line 52, before
        atoms.calc = XTB(method=args.method)
# after
        # Boundary condition comes from the structure, never from a default.
        assert all(atoms.get_pbc()[:2]), (
            f"{args.traj}: expected a,b periodic, got pbc={atoms.get_pbc()}")
        atoms.calc = TBLite(
            method=args.method,
            electronic_temperature=args.etemp,
            accuracy=1.0,
        )
```

and add the argument, so the electronic temperature is recorded rather than
inherited from whichever code is in use:

```python
    ap.add_argument("--etemp", type=float, default=300.0,
                    help="Electronic temperature in K, recorded in results JSON")
```

## 5. `inputs/frozen_optimiser_config.json`

A new version, since the boundary condition is part of the frozen method.
`triage_runs.py` hashes this file into every record, so bumping
`config_name` is what separates the old runs from the new ones.

```json
{
  "config_name": "irox_gfn2_bfgs_v4_periodic",
  "description": "Frozen optimiser settings, periodic. Supersedes irox_gfn2_bfgs_v3, which relaxed every slab as a finite cluster.",
  "xtb_version": "6.7.1",
  "tblite_version": "0.6.0",
  "engine": "tblite",
  "gfn_level": "GFN2-xTB",
  "pbc": [true, true, true],
  "vacuum_along_c_A": 30.0,
  "k_points": "gamma",
  "electronic_temperature_K": 300.0,
  "accuracy": 1.0,
  "ase_optimiser": "BFGS",
  "fmax_target_eV_per_A": 0.05,
  "max_steps": 200,
  "charge": 0,
  "multiplicity_handling": "xtb_default_uhf_0",
  "constraint": {
    "type": "FixAtoms",
    "rule": "freeze all atoms with z below z_freeze_A",
    "z_freeze_A": 21.01
  },
  "surface_threshold_A": 1.5,
  "random_seed": null
}
```

Two fields need a number from a test rather than from this file.
`vacuum_along_c_A` is a placeholder: xtb's GFN1 energy moved 0.718 eV when c
grew by 20 A, so run the vacuum convergence yourself with
`scripts/pbc_2d_vs_3d_check.py` under whichever engine you adopt and set the
converged value. `k_points` records gamma because that is what tblite 0.6.0
offers; if you need sampling, that is an argument for DFTB+ or CP2K.

## 6. One guard that prevents the whole class of bug

The slab builder wrote `(True, True, False)`, the state generator wrote
`(False, False, False)`, the triage gate assumed `(True, True, False)`, and
nothing compared them. Add a single check where structures enter the
optimiser, in `scripts/optimization_iro2_test.py` after the structure is
read:

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

with a matching `--config` argument defaulting to
`inputs/frozen_optimiser_config.json`.

## 7. Regenerate, do not patch

Existing relaxed structures were produced under the old boundary condition
and cannot be repaired by relabelling. Regenerate the states and rerun the
relaxations after changes 1 to 3 land. Do not overwrite `site_intended` at
any point, and keep the old outputs under the `irox_gfn2_bfgs_v3` config name
so the two sets stay distinguishable.

## New scripts added by this work

| File | Status |
| --- | --- |
| `scripts/pbc_feasibility_probe.py` | run on macOS arm64, results in `docs/pbc_feasibility_probe_macos_arm64.json` |
| `scripts/pbc_2d_vs_3d_check.py` | run on macOS arm64, results in `docs/pbc_2d_vs_3d_check_macos_arm64.json` |
| `scripts/singlepoint_edge_vs_interior.py` | NOT RUN |
| `scripts/relax_periodic.py` | NOT RUN |
