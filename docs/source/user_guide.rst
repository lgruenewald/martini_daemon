Installation
============


Running a simulation
====================

Simulations in Martini Daemon are ran using python run scripts. These
contain calls to Martini Daemon’s API, specifying simulation parameters
and input files. A simple energy minimization run script is provided
below, which should be adjustable to meet various needs.

::

   #!/usr/bin/env python3

   from martini_daemon.simulation import Simulation
   from martini_daemon.reporters.variables_reporter import VariablesReporter
   import openmm as mm

   # equilibration
   eq = Simulation(
     # input files
     top_path="system.top", geom_path="system.gro",
     # output files will be named after this
     sim_name="eq",
     # total MD steps to take
     md_steps=400000,
     # how often to write to the .xtc
     traj_frequency=5000,
     # no reactions during equilibration
     dm_frequency=0,
     # report thermodynamic variables, other reporters go here too...
     reporters=[
       # each reporter typically writes its own file with its own
       # specific extension
       VariablesReporter()
     ],
     # integrator+temperature coupling in one
     integrator=mm.LangevinMiddleIntegrator(
       # temperature
       298 * mm.unit.kelvin,
       # collision frequency
       1.0 / mm.unit.picosecond,
       # timestep
       0.02 * mm.unit.picosecond
     ),
     # coupling
     coupling=[
       mm.MonteCarloBarostat(
         # desired pressure
         1.0 * mm.unit.bar,
         # reference temperature
         298 * mm.unit.kelvin
       )
     ],
     # use CUDA with nvidia GPUs
     platform="CUDA"
     # Run on GPU 1 only
     context_parameters={"DeviceIndex": "1"}
   )
   # minimize energy first, saving the minimized coordinates to min.gro
   eq.minimize_energy(out="min.gro")
   # generate velocities at 298 K
   eq.generate_velocities(298)
   # run equilibration
   eq.simulate()

Selecting GPUs for the simulation can be done using context parameters,
as seen in the example above. Selecting CPU cores for the simulation can
be done with the ``taskset`` command. For example,
``taskset -c 0-31 ./run.py`` will limit run.py to cores 0 to 31.


Extracting the OpenMM system
============================

In case you want to use Martini Daemon as a Gromacs .top file parser to
run your (martini) simulations in OpenMM, you might want to just obtain
an OpenMM system, rather than use the abstraction layer provided on top
of it here. TODO -- write up

Including reactions
===================

A rough workflow for adding a reaction consists of several steps. First,
the desired reactions should be broken down to a mechanism, that can be
modelled. For each mechanistic step, the reactant and product molecules
should be parametrized in Martini. The difference between the two should
be written down as a list of new interactions, as well as old
interactions to break.

Second, a graph for the reactant should be constructed. Each graph match
will be added to a list of active known reactants (also called
fragments). Graphs specify a list of beads to match for, including their
names and types. It can be useful to name beads that react a special
way, which is then matched in the graph. Other atoms of importance, such
as those involved in angle or dihedral conditions, or for the formation
of new angles in the reaction should also be added to the graph. They
can be filtered using names, types as well. The interactions connecting
the beads in the graph should be specified. They can be specified in a
generic way, such as ``bond``, which matches every type of bond.
Alternatively, a specific filter can be used, such as ``harmonic_bond``
for bond type 1 only. If ambiguities arise, for example this often
happens with virtual site matches, bond type 5 ``connection`` can be
used to disambiguate graphs without changing the dynamics of the system.
Martini Daemon will keep an exhaustive and up to date list of all graph
matches within the system during the entire simulation. This fact should
be used for the removal of reacted reactants. If the reaction renames or
changes the type of the reacting bead, the graph will no longer match,
removing it from the list of potential reactants during a reaction.
Alternatively, if applicable, the graph should include a forbidden atom,
representing the other reactant that is bonded after the reaction.

Graphs start with the [frag] directive. In this directive, the first
word of each line specifies what information that line contains. First,
a unique name should be given to the graph, using the ``name`` keyword.
Second, a list of atoms should be given using the ``atom`` keyword. This
keyword has two alternate variations: ``atom?`` for optional beads and
``atom!`` for forbidden beads. After each of these directives, a name
has to be given to the atom (containing numbers or letters), followed by
the name and type filters.

::

   [frag]
   name example
   atom reactive BR1 SC6
   atom angle    C2  SC5
   bond reactive angle

The exhaustive graph matching has a consequence. If there is any type of
symmetry within the graph, all permutations with different system atom
to graph atom mapping will be separately added to the graph. This can be
undesirable, and can be prevented using the ``equivalent`` keyword. In
each equivalent line, all graph atoms specified will be considered
equivalent. When exchanged, equivalent atoms do not produce new graph
matches, removing the redundant matches.

The FragCountReporter will report the number of each fragment in the
system, for every frame. This can be useful to check, whether the
reactant fragments have the correct number of matches. It can be
imported as:

::

   from martini_daemon.reporters.topstar import FragCountReporter

Then, it has to be added to the list of reporters. This reporter also
adds a live number of fragments to the CLI during the simulation.

Third, the reaction template should be made. These are specified within
the ``[reaction]`` directive. Within the main body of this directive,
the reaction name must be specified.

::

   [reaction]
   example_reaction

After the reaction name, the list of reactants (graph names) should be
given. Reactions of only one fragment should specify one reactant here,
while reactions between two fragments should specify two.

::

   [reactants]
   example example

After this, the list of reaction conditions should be specified. A
reaction will only be accepted, if all specified conditions are true.
See the commented example below. Within this block, atoms are referenced
using a double index. The first part of each atom index specifies which
reactant the atom is a part of, while the second specifies the name
within the graph. For the possible conditions, see the example below, or
consult the reference.

::

   [conditions]
   ; Note: using the example graph specified above, which contains two atoms
   ; named reactive and angle respectively

   r_max 1:reactive 2:reactive 0.5 ; maximum distance of 0.5 between the main reactive beads
   r_min 1:reactive 2:reactive 0.3 ; minimum distance of 0.3
   angle_between 1:angle 1:reactive 2:reactive 70 100 ; an angle condition between 70 and 100 degrees
   angle_between 2:angle 2:reactive 1:reactive 70 100 ; same requirement, but from the other molecule's point of view
   dihedral_between 1:angle 1:reactive 2:reactive 2:angle 50 100 ; a dihedral condition, between 50 and 100 degrees

As a final part of reaction templates, the list of changes to topology
during a reaction should be specified. New interactions can be given
just like in moleculetypes. With the important difference, that atoms
need the double indexing here.

::

   [bonds]
   1:reactive 2:reactive 1 0.4 1000  ; a new harmonic bond formed in the reaction

   [angles]
   1:angle 1:reactive 2:reactive 1 85 100
   2:angle 2:reactive 1:reactive 1 85 100

Atom properties can be changed using the ``[redefine]`` directive. In
each line of this directive, a keyword specifies what is changed (name,
type, mass, charge), followed by the new value.

::

   [redefine]
   1:reactive name BRD type SN1
   2:reactive name BRD type SN1

Existing interactions can be removed too, using the ``[break]`` and
``[update]`` directives. Each line in these directives contains a list
of double indexed atoms, in which the interactions should be broken
according to the rules of the directives.

::

   ; for a different reaction containing two atoms, left and right
   [break]
   1:left 1:right

The break directive completely disconnects the molecule between two
specified atoms. The bonds, angles, dihedrals, exclusions that contain
both atoms will all be removed. The update directive has the same form
as the break directive, but is more careful. It only breaks interactions
that are fully between the specified atoms. So, in the case of two atoms
specified, it will break the bond and exclusion, but not angles or
dihedrals that contain them and other atoms as well.

Lastly, the setup should be tested. It’s recommended to run simulations
with the ReactionReporter, as it reports each reaction, along with the
frame it happens, the internal fragment ID of reactants, and a list of
atoms within the fragment. It can be imported as:

.. code:: py

   from martini_daemon.reporters.topstar import ReactionReporter

Then, it has to be constructed and added to the list of reporters within
the simulation. This reporter also adds a cumulative reaction counter to
the CLI during the simulation.

Specific simulation requirements
================================

It is possible to temporarily reduce the timestep of the simulation
after reactions, if required using “reaction sensitive” integrators. See
example for below.

::

   #!/usr/bin/env python3

   from martini_daemon import simulation
   from martini_daemon.reporters.bond_reporter import BondReporter
   from martini_daemon.reporters.topstar import ReactionReporter
   from martini_daemon.components.reaction_sensitive_integrator import ReactionSensitiveLangevinIntegrator

   sim = simulation.Simulation(
       top_path="system.top", gro_path="system.gro",
       sim_name="out",
       reporters=[
           BondReporter(),
           ReactionReporter(molid=True),
       ],
       md_steps=100000000, dm_frequency=100,
       xtc_frequency=5000,
       # 0.02 ps timestep, 298 K, 1 ps^-1 friction
       # subdivision of 4, for 5 timesteps
       integrator=ReactionSensitiveLangevinIntegrator(0.02, 298, 1., 4, 5),
   )
   sim.minimize_energy()
   sim.generate_velocities(300)
   sim.simulate()

The subdivision has to be a positive integer. During the short
equilibration, each timestep is divided into this many sub-time steps.
In this example, this means 5 fs timesteps. The number of steps in this
example is 5, which get divided into 20 fs timesteps. It is ensured,
that XTC frames remain evenly spaced, regradless whether there are
reactions happening.

Using this may have a performance impact, so it is recommended to try to
specify reaction conditions that do not require this short post-reaction
equilibration.

Reporting
=========

There are various other reporters which can be useful worth mentioning
briefly in this guide.

Variables Reporter
------------------

Reports thermodynamic variables, such as kinetic, potential and total
energies, temperature and box size. A new entry is written every time
the XTC trajectory file is written, for easy analysis.

.. code:: py

   from martini_daemon.reporters.variables_reporter import VariablesReporter

Checkpoint Reporter
-------------------

The checkpoint reporter is constructed with a specific interval, at
which it will write simulation checkpoints. A simulation can be
continued from these checkpoints, when using the exact same martini
daemon version.

.. code:: py

   from martini_daemon.reporters.checkpoint_reporter import CheckpointReporter

Checkpoints can be loaded by specifying the chk_path argument of
Simulation instead of specifying the geom_path and top_path arguments.

Atom Reporter
-------------

The atom reporter writes atom information (name, type, charge, mass) for
each XTC frame.

.. code:: py

   from martini_daemon.reporters.atom_reporter import AtomReporter, read_atoms

   ...

   n_frames, natoms, names, types, charges, masses = read_atoms("out.atoms")

Bond Reporter
-------------

The bond reporter writes a list of bonds in the system for each XTC
frame. All entries in the ``[bonds]`` directive, constraints are
considered bonds. Virtual sites are considered bonds too, between the
virtual particle and all constructing particles respectively.

.. code:: py

   from martini_daemon.reporters.bond_reporter import BondReporter, read_bonds

   ...

   n_frames, n_atoms, bond_frames = read_bonds("out.bonds")

Helpers
=======

The helpers folder contains helpers that facilitate analysis or
visualization.

VMD
---

The VMD helper, together with the VMD script located in the ``tcl``
directory of this repository facilitate the visualization of bonds
during trajectories with bond formation and breakage. The workflow is as
follows:

1. Have a ``.bonds`` BondReporter output and an XTC file, which is
   desired to be visualized. If molecules are made whole across the PBC,
   or similar transformations, those should be done first.
2. Generate the file that is read by the VMD script from the ``.bonds``
   and ``.xtc`` files.

.. code:: py

   from martini_daemon.helpers.vmd import generate_vmd_readable_bonds
   from martini_daemon.bond_reporter import read_bonds

   n_frames, n_atoms, bond_frames = read_bonds("out.bonds")
   generate_vmd_readable_bonds(
     bond_frames,
     "out.xtc",  # XTC path
     "out.z"  # output path
   )

3. (Optional, for large trajectories) look at the README in the ``tcl``
   directory of this repo for instructions to compile the shared library
   helper.
4. Make a visualization script that loads the tcl script and load the
   bond trajectory.

.. code:: tcl

   # (optional) if step 3 was completed:
   load /path/to/martini_daemon/tcl/bond_loader.so
   # always mandatory:
   source /path/to/martini_daemon/tcl/daemon.tcl

   # loads the .gro and .xtc file, deleting the extra frame from the .gro
   daemon_open out.gro out.xtc
   # loads the bond list generated in generate_vmd_readable_bonds
   daemon_bonds out.z

   # ... other commands to set visualization state

This script can then be loaded using the ``-e`` flag of VMD. If the
script is saved as ``vis.tcl``, it can be opened as:

::

   vmd -e vis.tcl
