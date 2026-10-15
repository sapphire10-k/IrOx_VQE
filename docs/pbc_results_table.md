# Results table: step 1 feasibility

Unrelaxed canonical clean slab, Ir32O64, 96 atoms. Cell 12.72311 x 6.30880 x 42.01316 A, all angles 90 deg. Slab thickness 12.0132 A, vacuum along c 30.0000 A.

Environment: macOS-26.6-arm64-arm-64bit, Python 3.11.17, ASE 3.29.0, NumPy 2.4.6, xtb 6.7.1, xtb-python 22.1, tblite 0.6.0.

## Table 1. Single point, by engine, method and boundary condition

| Engine | Method | pbc | Status | Energy (eV) | fmax all atoms (eV/A) | Wall (s) |
| --- | --- | --- | --- | --- | --- | --- |
| xtb-python | GFN2-xTB | FFF | ok | -10376.841415 | 4.0802 | 3.824 |
| xtb-python | GFN2-xTB | TTF | failed (multipoles not available with PBC) | - | - | 0.149 |
| xtb-python | GFN2-xTB | TTT | failed (multipoles not available with PBC) | - | - | 0.152 |
| tblite | GFN2-xTB | FFF | crashed (SIGSEGV) | - | - | - |
| tblite | GFN2-xTB | TTF | crashed (SIGSEGV) | - | - | - |
| tblite | GFN2-xTB | TTT | crashed (SIGSEGV) | - | - | - |
| xtb-python | GFN1-xTB | FFF | ok | -11618.104779 | 6.1859 | 2.42 |
| xtb-python | GFN1-xTB | TTF | ok | -11693.525683 | 1.2063 | 2.939 |
| xtb-python | GFN1-xTB | TTT | ok | -11693.525683 | 1.2063 | 3.071 |
| tblite | GFN1-xTB | FFF | ok | -11618.104787 | 8.8366 | 2.444 |
| tblite | GFN1-xTB | TTF | ok | -11683.433681 | 1.3282 | 2.579 |
| tblite | GFN1-xTB | TTT | ok | -11683.433681 | 1.3282 | 2.678 |

## Table 2. Vacuum sensitivity, GFN1-xTB, atoms held fixed in Cartesian space

| Engine | pbc | c = 42.013 A (eV) | c = 62.013 A (eV) | dE (eV) | Reading |
| --- | --- | --- | --- | --- | --- |
| xtb-python | TTF | -11693.525682659 | -11694.243624371 | -0.717941712 | c is used: this is 3D periodicity |
| xtb-python | TTT | -11693.525682659 | -11694.243624371 | -0.717941712 | c is used: this is 3D periodicity |
| tblite | TTF | -11683.433680727 | -11683.433681263 | -0.000000536 | c is irrelevant: consistent with true 2D |
| tblite | TTT | -11683.433680727 | -11683.433681263 | -0.000000536 | c is irrelevant: consistent with true 2D |

## Table 3. Cross-checks derived from Table 1

| Quantity | Value |
| --- | --- |
| GFN1 non-periodic, xtb-python minus tblite | +8.30e-06 eV |
| GFN1 periodic (TTT), xtb-python minus tblite | -10.0920 eV |
| GFN1 open to periodic, xtb-python | -75.4209 eV |
| GFN1 open to periodic, tblite | -65.3289 eV |
| GFN1 fmax, open to periodic, xtb-python | 6.1859 to 1.2063 eV/A |
| GFN1 fmax, open to periodic, tblite | 8.8366 to 1.3282 eV/A |

## Steps 2 to 4

NOT RUN. `scripts/singlepoint_edge_vs_interior.py` (step 2) and `scripts/relax_periodic.py` (steps 3 and 4) are written and syntax checked but unexecuted. Steps 3 and 4 as originally specified, a periodic GFN2 relaxation, cannot run until the tblite question is settled on Linux.

