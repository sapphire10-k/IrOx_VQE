#Periodic boundary conditions for H on Iridium Oxide (110): recommendation
For discussion with Prof David Benoit. Prepared 14 October 2026.

##The decisive result
GFN2-xTB cannot be run with periodic boundary conditions in xtb 6.7.1. The library refuses, with its own message:


-2- xtb_calculator_singlepoint: Electronic structure method terminated
-1- scf: Multipoles not available with PBC
This is a property of the Hamiltonian as implemented, not a missing flag. GFN2-xTB carries atomic dipoles and quadrupoles, and that anisotropic electrostatic expansion has no periodic implementation in xtb. GFN1-xTB, whose electrostatics is monopole only, runs periodically without complaint. The refusal is identical for (True, True, False) and (True, True, True).

A second route, tblite 0.6.0, could not be tested here: its GFN2 path segfaults on this machine (macOS arm64) even for a three atom water molecule with no periodicity at all, while its GFN1 path works normally. That is a build fault on this platform and says nothing about periodicity. It is the one open question, and it should be settled on your Linux machine by running scripts/pbc_feasibility_probe.py there.

Two further findings that change how the question should be asked
pbc = (True, True, False) is not honoured as two dimensional periodicity by either engine. Both return bit identical energies for (T,T,F) and (T,T,T). Growing the c axis by 20 A with the atoms untouched moves the xtb energy by 0.718 eV under both flags, so xtb is using the c vector and (T,T,F) is three dimensional periodicity in disguise. tblite's energy is flat to 5e-7 eV under the same change, so tblite ignores c entirely. There is therefore no 2D option to choose: you supply a 3D cell with vacuum and accept what the code does with the third direction. The 0.718 eV drift also says the canonical 30.000 A vacuum is not converged for xtb and would need testing.

The two engines disagree by 10.09 eV on the same periodic GFN1 calculation, having agreed to 8.3e-6 eV on the same structure without periodicity. The non-periodic agreement is close enough to show that their defaults, electronic temperature included, already match, so the 10 eV gap lies in the periodic treatment itself. Periodic xTB is not yet an interchangeable, settled method between codes, and whichever is used must be named and versioned in the methods section.

##Does periodicity remove the artefact?
The evidence supports your hypothesis, with the caveat that it is GFN1 and not GFN2. On the unrelaxed canonical clean slab the maximum force falls from 6.19 eV/A (xtb) and 8.84 eV/A (tblite) with open boundaries to 1.21 and 1.33 eV/A with periodicity, and the total energy falls by 65 to 75 eV. Forces of 6 to 9 eV/A at step zero mean a BFGS trajectory dominated by relieving missing neighbours at the free faces. The mechanism is geometric rather than Hamiltonian specific, so it will apply to GFN2 as well.

On keeping non-periodic relaxations with a cluster cut
I do not think this is defensible for this cell, and the objection is geometric rather than a matter of degree.

The b axis is 6.3088 A and the a axis is 12.7231 A. With the boundaries open, all four a and b faces are free surfaces, so the deepest point in the cell, its centre, is only min(a,b)/2 = 3.1544 A from the nearest face. No atom in this slab can be more than 3.1544 A from a free surface. Against your own 2.30 A Ir-O neighbour criterion, excluding even one coordination shell around each face leaves a region under 1 A thick. There is no artefact free interior to cut to. Your audit says the same thing from the other direction: 36 of 64 O atoms lose Ir neighbours, which is 56% of the oxygen sublattice, not a rim.

Three of the four adsorption sites make this concrete. Reading the canonical input, o21 sits at (0.0000, 0.0000), o45 at (0.0000, 3.1544) and o69 at (6.3616, 0.0000), all on a cell boundary plane. Only o93, at (6.3616, 3.1544), is at the cell centre in both a and b. A cut that excluded edge artefacts would exclude three of the four sites under study.

The deeper objection is that the cut is applied after relaxation. The final geometry everywhere, o93 included, was produced by an optimiser responding to forces that included fictitious edge chemistry. Deleting the visibly broken atoms afterwards does not undo their influence on where the survivors came to rest. Your H*-o69 run shows exactly this: its short O-O contacts are 5.05 and 9.22 A from the site O, so the site looked locally healthy while the cell around it reconstructed.

To make the cut defensible you would have to demonstrate that the site local geometry is converged with respect to distance to the nearest free face, by relaxing progressively larger laterally expanded cells and showing the geometry within R A of the site stops moving. Allowing two intact coordination shells between a 5 A cluster and any free face implies lateral dimensions of order 26 to 30 A, roughly a 4x4 expansion and about 1500 atoms. That is more expensive than the periodic calculation it was meant to avoid.

##Recommendation
Run the feasibility probe on Linux first. It is seconds of compute and it decides everything below. If tblite gives periodic GFN2 there, take it: the Hamiltonian, the frozen optimiser settings and the whole downstream pipeline stay as they are, and only the calculator and the pbc flag change.
If periodic GFN2 is unavailable, move to periodic DFT. The repository is already part way there: inputs/pseudopotentials, scripts/optimization_iro2_qe_pseudopotential.py, and a canonical slab that is already a Quantum ESPRESSO input with CELL_PARAMETERS. For a 96 atom Iridium Oxide slab this is the route least likely to attract a referee's objection. It will need more than gamma point sampling.
Treat periodic GFN1-xTB as a screening tool only, not as the production geometry source. It works today, but it changes the Hamiltonian and carries the unresolved 10 eV inter-code disagreement.
Do not report geometry from the current open boundary relaxations, H*-o69 included, while this is unresolved.
What each choice means for the methods section
Periodic GFN2 via tblite. The cleanest story. You state the code and version, the 3D cell with its vacuum thickness, gamma point sampling and a vacuum convergence test. You must also note that H* is now an ordered overlayer, with adsorbate images 6.3088 A apart along b, rather than an isolated adsorbate. That is a change in what is being computed and a referee will look for it.

Periodic DFT. The conventional choice, and the methods section writes itself: functional, pseudopotentials, cutoffs, k-point mesh, smearing, dipole correction for the asymmetric H covered slab. The cost moves to the cluster.

Periodic GFN1. Requires you to justify dropping to monopole electrostatics for a transition metal oxide surface, name the implementation, and show the vacuum convergence. I would expect questions.

Non-periodic with a cluster cut. Requires the lateral convergence study described above. Without it the geometry is not defensible, and with it the cost exceeds the periodic route.

What was and was not measured
Produced here, on the unrelaxed canonical clean slab, macOS arm64, xtb 6.7.1, xtb-python 22.1, tblite 0.6.0, ASE 3.29.0, NumPy 2.4.6, Python 3.11.17: the feasibility matrix, the vacuum sensitivity test, and the geometric facts read from the canonical input.

Not run: steps 2, 3 and 4. scripts/singlepoint_edge_vs_interior.py and scripts/relax_periodic.py are written and reviewed but unexecuted, and both carry an EXECUTION STATUS header saying so. Note that steps 3 and 4 as originally specified, a periodic GFN2 relaxation, cannot run until the tblite question is settled.
