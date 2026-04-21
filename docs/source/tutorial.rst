Tutorial: Thiol disulfide bond formation
========================================

Authors: Aster Kovacs, Linus Grunewald


Martini Daemon is a tool to perform template based chemical reactions in molecular dynamics simulations with the
Martini force field. This page contains a step-by-step tutorial that walks the user through
the intended Martini Daemon workflow step by step.

This tutorial focuses on reproducing the Reactive Martini model (https://pubs.acs.org/doi/10.1021/acs.jctc.2c01186)
within Martini Daemon. Reactive Martini is based on tabulated potentials and GROMACS, to include a reaction template
which allows for the formation of disulfide bridges from thiols. Benzene-1,3-dithiol monomer molecules are simulated to
form macrocycles of sizes mainly 3 and 4, that is cyclic trimers and tetramers.

In this tutorial, a similar, but more advanced system is constructed. The reaction template introduced in this
tutorial will allow a thiolate anion to cleave a pre-existing disulfide bridge, displacing a thiolate anion and forming
a new disulfide bridge. This mechanism would be previously difficult and impractical to precisely model with
tabulated potentials. Martini Daemon allows for the precise definition of the desired reaction mechanism.
Do note, that the reaction mechanism becomes another *input* to the simulation. with Martini Daemon, you can only study
the impact of the reaction on larger systems, rather than the individual elementary reaction steps themselves.

Pre-requisites
--------------

The first step in the Martini Daemon workflow is the creation of Martini models for the reactant and product
molecules. Help for this step can be found on the https://cgmartini.nl/tutorials website, especially the small
molecule parametrization.

An important consideration is to try to map the reactants and products to a similar
bead mapping. Martini Daemon cannot change the number of beads during a reaction, and the output file formats (.xtc),
or analysis libraries (mdtraj, MDAnalysis) are also not well suited for simulations with a variable number of atoms.
Changing atom types, charges, mass, forming and breaking bonds, angles, dihedrals and adding/removing exclusions is
supported. The difference between the reactant(s) and product(s) should be implementable with only those changes.

For this tutorial, 3,5-disulfanylbenzoic acid will be used as the monomer. Here is the ``.itp`` file for the monomer.

::

    [ moleculetype ]
    DSBO 1

    [ atoms ]
     1 SQ3  1 DSBO S1   1 -1.0 66.0
     2 SQ3  1 DSBO S2   2 -1.0 66.0
     3 TC5  1 DSBO BNZ1 3  0.0  0.0
     4 SQ5n 1 DSBO CAC1 4 -1.0 66.0
     5 U    1 DSBO VS1  5  0.0  0.0

    [ constraints ]
     1  2   1  0.38882
     1  4   1  0.48533
     2  4   1  0.48549

    [ virtual_sites3 ]
     3 1 2 4 1  0.25552  0.48859

    [ virtual_sitesn ]
     5 1 1 2 3

    [ exclusions ]
     1 2 3 4
     2 3 4
     3 4

The dimer was parametrized as:

::

    [ moleculetype ]
    DSB2 1

    [ atoms ]
     1 SQ3  1 DSB2 S1   1  -1.0 72.0
     2 SC6  1 DSB2 S2   2  0.0  72.0
     3 TC5  1 DSB2 BNZ1 3  0.0  0.0
     4 SQ5n 1 DSB2 CAC1 4  -1.0 72.0
     5 U    1 DSB2 VS1  5  0.0  0.0
     6 SC6  1 DSB2 S3   6  0.0  72.0
     7 SQ3  1 DSB2 S4   7  -1.0 72.0
     8 TC5  1 DSB2 BNZ2 8  0.0  0.0
     9 SQ5n 1 DSB2 CAC2 9  -1.0 72.0
    10 U    1 DSB2 VS2  10 0.0  0.0

    [ bonds ]
     6  2   2  0.294 18000; S3_S2

    [ constraints ]
     1  2   1  0.38882
     1  4   1  0.48533
     2  4   1  0.48549
     6  7   1  0.38882
     6  9   1  0.48533
     7  9   1  0.48549

    [ angles ]
     5  2 6 10 120 220 ; VS1_S2_S3
     10 6 2 10 120 220 ; VS2_S3_S2
     1  5 6 10 110 5 ; S1_VS1_S3
     7 10 2 10 110 5 ; S4_VS2_S2

    [ dihedrals ]
     5 2 6 10 2 -70 90 ; VS1_S2_S3_VS2
     3 1 5 6 3  -0.002 0.011 0.024 -0.065 -0.019 0.053 ; BNZ1_S1_VS1_S3
     8 7 10 2 3  -0.002 0.011 0.024 -0.065 -0.019 0.053  ; BNZ2_S4_VS2_S2

    [ virtual_sites3 ]
     3 1 2 4  1  0.25552  0.48859
     8 6 7 9  1  0.25552  0.48859

    [ virtual_sitesn ]
      5 1 1 2 3
     10 1 6 7 8

    [ exclusions ]
    1 2 3 4
    2 3 4
    3 4
    6 7 8 9
    7 8 9
    8 9

The additional bond, angles and dihedrals connecting the two monomer units are indicated with comments.

Defining reactants
------------------

Defining reaction templates
---------------------------

Running the simulation
----------------------

Analyzing the output
--------------------
