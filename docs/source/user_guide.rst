User Guide
==========

This page of the documentation gives a quick primer on various topics.

Introduction to run scripts
---------------------------

Simulations in Martini Daemon are ran using python run scripts. These
contain calls to Martini Daemon’s API, specifying simulation parameters
and input files. A simple equilibration run script is provided
below, which should be adjustable to meet various needs.

::

   #!/usr/bin/env python3

   from martini_daemon import Simulation, VariablesReporter, XTCReporter
   import openmm as mm

   eq = Simulation(
     # path to Gromacs Topology
     top_path="system.top",
     # path to Starting geometry
     geom_path="system.gro",
     # output filenames will be prefixed with this name
     sim_name="eq",
     # total MD steps to take
     md_steps=400000,
     # how often to write to the .xtc
     traj_frequency=5000,
     # no reactions during equilibration
     dm_frequency=0,
     # report thermodynamic variables, other reporters go here too...
     reporters=[
       # each reporter typically writes its own file with its own specific extension
       VariablesReporter(),
       XTCReporter(),
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
       ),
       mm.CMMotionRemover()
     ],
     # It's recommended to use CUDA with nvidia GPUs for optimal performance
     # see the available platforms by running python3 -m openmm.testInstallation
     platform="CUDA"
   )
   # minimize energy first, saving the minimized coordinates to min.gro
   eq.minimize_energy(out="min.gro")
   # generate velocities at 298 K
   eq.generate_velocities(298)
   # run equilibration
   eq.simulate()

See :doc:`/autoapi/martini_daemon/Simulation` for the whole list of available arguments, attributes and methods.

Including reactions - TODO rewrite
-------------------

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

Reporters
---------

To produce output from simulations, reporters have to be added to simulations.
Reporters are Python class instances that inherit the ``Reporter`` base class,
and their methods get called by ``Simulation`` at specified events during
simulations.

Reporters can be broadly divided into two categories:

* Some perform reporting at a pre-defined trajectory frequency:
    * This frequency is the ``traj_frequency`` argument to the constructor of ``Simulation``.
    * This frequency is the same per-simulation to make analysis easier -- frame indices are synchronized across reporters.
    * For a successful simulation, there is a frame written at the start, when the current MD step modulo frequency is 0, and at the very end of the simulation.
    * :doc:`/autoapi/martini_daemon/VariablesReporter` reports thermodynamic variables for each trajectory frame.
    * :doc:`/autoapi/martini_daemon/XTCReporter` reports atom positions and the periodic box for each trajectory frame.
    * :doc:`/autoapi/martini_daemon/ToptrajReporter` creates a "topology trajectory", reporting atom properties and bonds for each trajectory frame.
    * :doc:`/autoapi/martini_daemon/FragCountReporter` reports the number of fragments at each trajectory frame. It also adds the current total number of fragments to the interactive line on stdout.

* Some perform reporting related to reactions happening in the system:
    * :doc:`/autoapi/martini_daemon/LocalMinimizer` is not a traditional reporter. It locally minimizes the energy after reactions.
    * :doc:`/autoapi/martini_daemon/ReactionReporter` logs all reactions and reactants to a file as they happen.
    * :doc:`/autoapi/martini_daemon/ReactionEnergyReporter` reports thermodynamic variables before and after reactions. Optionally, it can write coordinates too, which can be helpful to debug local minimizations.

Of course, technically a reporter can do something at both, it is up to the implementation to choose which
callbacks to hook on. User-defined reporters are supported, they can use the public API of simulation
and all its public attributes to perform tasks during callbacks. Callbacks are called in the order
they are passed to Simulation. See :doc:`/extending` for more detail.

Analysis - TODO
--------

Visualization
-------------

The ``.toptraj`` files generated by Martini Daemon's :doc:`/autoapi/martini_daemon/ToptrajReporter` can be
used to dynamically visualize the topology as it changes during a trajectory.

Currently, there is a plugin for `VMD`_, which can be found in the ``/vmd_plugin`` folder of the Martini Daemon
git repository. This plugin is written in C, so it needs to be compiled first.
Instructions on how to build and use this plugin are in the README in said folder, as it can be viewed as a separate
component.
