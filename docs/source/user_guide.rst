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
   # write final coordinates, velocities and box to a .gro file
   eq.save_geometry("eq.gro")

See :doc:`/autoapi/martini_daemon/Simulation` for the whole list of available arguments, attributes and methods.

Including reactions
-------------------

A rough workflow for adding a reaction consists of several steps.

1. The desired reactions should be broken down to a mechanism.
2. For each mechanistic step, the reactant and product molecules should be parametrized in Martini.
3. The difference between the two should be written down as a list of new interactions, as well as old interactions to break.
4. A graph for the reactant should be constructed.
5. A reaction template should be defined.

Creating reactant graphs
------------------------

The graph matching algorithm provides a powerful method to create
a "contract" which defines a reactant in a flexible manner.
Each graph match will be added to a list of active known reactants (also called
**fragments**). Graphs specify a list of beads to match for, including their
names and types. Graphs also specify interactions to match for, which connect
the beads. These interactions can be generic, such as ``bond``, or specific,
such as ``harmonic_bond`` or ``morse_bond``. All graph beads should be connected
to each-other with such interactions.
The basic syntax for defining graphs is done using the ``[graph]`` directive. Each line within this directive
starts with a keyword, and then parameters that follow. All graphs must be given a name, using the ``name`` keyword.
Atoms should be defined with the ``atom`` keyword, followed by the atom name, the atom name filter, and the atom
type filter.

Let us consider a simple example,
where a bead representing an alcohol, called ROH with type SP3 is matched. Another bead is included.
This bead will be used for angle conditions and angle forces. Its graph node name will be tail,
with a name match pattern that matches any atom name starting with C and any type.
A generic bond interaction is matched, so they can be connected using any of the available bond types without
needing to update the graph.

::

    [graph]
    name alc
    atom oh   ROH SP3
    atom tail C* *
    bond oh tail

If this graph is specified, Martini Daemon will find all beads with the name ROH, type SP3, that are bonded to any
bead whose name starts with C, and it will collect all the matches in the fragment list.
The graph matching algorithm will keep an always up to date list of all graph matches for the detection algorithm.
This means, that once the **contract** specified by the graph description above is broken, the corresponding
fragment is removed. This is done by matching all atoms at the start of the simulation,
and then re-calculating matches for atoms who participate in reactions, as well as their "neighbors", under
the assumption that only reactions change the topology of the system, and only topology changes can invalidate
graphs.

Interaction filters with more than two atoms are also supported. Angles or dihedrals can be matched as well.
Interaction filters such as ``angle oh bead2 bead3`` are valid. This specific one will match any angle consisting
of the three specified beads, in any order.

Important! Graph matches are considered valid only if all atoms and interactions specified are matched. If a single
atom or interaction is missing, the graph match is considered invalid. Every interaction can only be matched once.
This can lead to complications when making complex graphs that try to match multiple interactions that include
more than two atoms, all with the same filter. While in some cases they may be technically possible,
such graphs are advised against, and are not supported.

In All-Atom models, reactive functional groups can be easily identified based on elements only
(-OH, -NH2, -COOH). This lends itself to well-defined reactant graph definitions without further complexity.
With Coarse Grained models, such as Martini, entire functional groups often correspond to a single specific bead.
The information of the true identity of the functional group remains with whoever is familiar with the model,
rather than it being implicit in the atom types alone. This is why atom names are matched as well, rather than
atom types only. Dynamic matching with ``*`` is still allowed to retain some dynamic behavior, as atom names
often have various suffixes. ``?`` is also allowed, which matches a single character. Curly braces allow
matching one character from a subset of characters. ``{123}`` matches a single character that can be ``1``, ``2``
or ``3``. It remains up the user to specify graphs that are general enough, for ease of use, but strict enough,
so they match only what was intended.

An important property of the graph matching algorithm, that it tries to be exhaustive. It will consider all
possible matches, and add them all to the list of fragments. Two matches are considered identical if for each
matched graph atom, the same atom index within the system is matched. For our example, this means, that if there
are two different beads bonded to the alcohol bead, with names that start with C, two different matches can be obtained.
Both matches, by default, will be added to the fragment list as viable reactants that fulfill the contract specified
in the graph description. If this is undesirable, it is important to make graphs that do not have such ambiguity.
There are multiple strategies for this:

* Use stricter name and type filters than in this example. If possible, change the names of the atoms in the .itp
  if they represent fundamentally different beads.
* Include other atoms in the graph that can be used to disambiguate.
* Add connections (bond type 5) to the molecule, and match them using the ``connection`` filter.
  This type of bond does not add any force to it, it is used only as a unit of topological information.
  Of course, this also means, that visualizations might draw a bond between the two particles, which may be
  undesirable.
* If the two beads are fully identical, and it is truly arbitrary which one is which, it may be okay to allow
  both matches.

The consequence of this exhaustiveness also comes up when matching multiple identical beads at the same time.
Consider this modified example:

::

    [graph]
    name alc
    atom oh    ROH SP3
    atom tail  C*  *
    atom tail2 C*  *
    bond oh tail
    bond oh tail2
    equivalent tail tail2

The keyword ``equivalent`` is used to tell the graph matching algorithm, that the beads tail and tail2 are,
for all purposes, equivalent. The only change this keyword introduces is to the graph matching algorithm itself.
On exchanging the two atom indices, the algorithm will not consider it a new match. Without the equivalent keyword,
two matches would be found, with tail and tail2 exchanged, and both would be separately added to the fragment list.
Auto-detection of equivalence is currently not possible, because when the graphs are defined,
it is unknown how reaction templates use them, so this must be specified and verified manually.

Removing reactant graphs on reactions
-------------------------------------

It is important to define reactant graphs and reactions in a manner, such that during a reaction, they no longer will
be valid reactants that can still react in the same reaction (unless this is explicitly desired). This is done
implicitly in Martini Daemon, by specifying a graph that no longer matches after the reaction.

Here is a list of general strategies for achieving this:

* In reactions that remove bonds, this can be automatically achieved by specifying a required bond between the two
  bonded atoms.
* In reactions that change atom names or types, this can also automatically be achieved with atom name and type filters,
  if they no longer match after the reaction.
* In other cases, a powerful tool can be the inclusion of forbidden atoms in graph specifications. The rest of this
  subsection will give an introduction to them.

Forbidden atoms are specified similarly to regular atoms, but they start with the `atom!` keyword. These atoms
function identically to regular atoms, with one difference. With only regular atoms, a graph match is complete if
all atoms and interactions were matched. If a forbidden atom is included, the interactions to the forbidden atom
only specify how to find it. If such a forbidden atom is matched, the graph match is considered invalid. If such a
forbidden atom could not be matched, the graph match remains valid. This is enabled by the eagerness of the
graph matching algorithm. The graph matching algorithm only considers whether a partial match is valid or not, once no
more atoms can be added to it based on the filters.

Note, that there is no way to group forbidden atoms together -- that is, rules such as "it is forbidden to have a
bead of this name and type, that is bonded to this other bead with this name and type" are not possible. Another
way to think about this restriction is that forbidden atoms can only represent a "one deep" layer around the normal
atoms. This is important, as this means, that a graph can only become forbidden if a change occurred to its direct
neighbors. This means, that graphs only need to be recalculated if their direct neighbors change. In order to
prevent users from attempting to group forbidden atoms, there is an error message raised if there is an interaction
filter connecting a forbidden atom to another forbidden atom.

Note, forbidden atoms can show up as a ``-1`` in some reporter outputs, since they are represented with a -1
internally in the fragment list.

Optional atoms
--------------

Optional atoms allow for some extra complexity in some cases. They function similarly to forbidden atoms,
with the difference that the graph is valid regardless of whether they are there or not.

General rules for optional atoms:

* Optional atoms are defined with the ``atom?`` keyword, followed by graph node name, name filter, type filter.
* Missing optional atoms show up as ``-1`` in some reporter outputs, as they are represented with a -1 internally.
* If a reaction condition references a missing optional atom, the condition is ignored.
* If a modification template entry references a missing optional atom, the whole entry is ignored.
  An exception to this is ``[update]``, as the semantics of that specifically make more sense that way.
* Due to the eagerness of the graph matching algorithm, if an optional atom can be matched, only the graph match
  that includes the optional atom is added to the fragment list.
* Optional atoms cannot be "grouped", that is interaction filters within the graph can only contain at most one
  optional atom per interaction filter. They also cannot be grouped with forbidden atoms. This facilitates
  graph recalculation on direct neighbor change only.
* Due to neighbor recalculation, if an existing fragment with a missing optional atom suddenly has an optional atom
  available to it, it will be recalculated to include it.

Debugging graphs
----------------

The reporter :doc:`/autoapi/martini_daemon/FragmentsDump` is specifically designed to facilitate debugging
the graph matching algorithm, by printing the list of matched fragments every time there was any recalculation
of graph matches.

Reaction templates
------------------

Reaction templates define:

* Reactants that should react,
* Reaction conditions for the reaction algorithm,
* Topology modifications to perform.

All reaction templates should start with the ``[reaction]`` directive, containing the reaction name.

::

   [reaction]
   example_reaction

The list of reactants should be given as a list of graph names in the ``[reactants]``directive.
The fragments in the fragment list with this name will be considered during the detection algorithm.
Up to three reactants are supported.

::

   [reactants]
   alc alc

Afterwards, reaction conditions should be specified.

In reaction conditions and topology modifications, individual atoms of the reactant graphs can be referenced.
This is done by referencing the (1 based) index of the reactant, followed by the graph node name within that graph.
This indexing has to happen in a single token, so they are separated by a colon and no whitespace.

Some example conditions are given below, with comments as explanations. Each of the listed conditions must be
fulfilled in order to accept a reaction. Consult the :doc:`/reference` for a full explanation of possible keywords.

::

   [conditions]
   r_max 1:oh 2:oh 0.5 ; maximum distance of 0.5 between the beads named "oh" in the graph
   r_min 1:oh 2:oh 0.3 ; minimum distance of 0.3
   angle 1:tail 1:oh 2:oh 70 to 100 ; an angle condition between 70 and 100 degrees
   angle 2:tail 2:oh 1:oh 70 to 100 ; same requirement, but from the other molecule's point of view
   dihedral 1:tail 1:oh 2:oh 2:tail 50 to 100 ; a dihedral condition, between 50 and 100 degrees
   p 0.5 ; if all others are fulfilled, 50% probability of accepting it

The following aspects should be considered when making reaction conditions:

* If forming a bond, generally there should be a maximum distance requirement for the given pair of atoms.
    * Forming a bond far away from the equilibrium bond length may result in a too high potential energy right after
      the reaction, or may result in other molecules still being "sandwiched" between the two reactants.
* If removing a bond, it can be useful to have a minimum distance requirement for the given pair of atoms.
* If forming an angle, sometimes it is necessary to have an angle condition.
* Sometimes, to get the right product geometry, distance, angle or dihedral reaction conditions are needed.
    * Some multi-functional molecules may self-react in undesirable ways if such conditions are too broad.
* Geometry conditions reduce reaction rates. Generally angle conditions have a larger impact than distance conditions,
  while dihedral conditions have an even larger impact.
* When multiple reactants are specified, a maximum distance
  must be present to "connect" them. This is because the detection algorithm constructs a KdTree of all fragment
  locations in space, which is used to accelerate it. Without a maximum distance condition, all combinations
  of reactants would need to be considered. This distance can be the same distance used for forming the bond.
* When tuning the reaction rates, be aware that:
    * Not only reaction rates relative to other reactions have an effect, but also relative to all processes happening
      in the simulation, including diffusion.
    * Different geometry conditions will result in different reaction rates. It may be useful to set reaction conditions
      at geometries that are somewhat higher in the potential energy surface, emulating a sort of "activation energy".
    * If using the probability reaction condition, the same reaction probability may result in different rates based on
      the other conditions, as the probability is only applied "on top" of any other conditions.

Finally, the list of changes to topology during a reaction should be specified. Currently, this is a quite manual
process. Users should first write down the difference between the reactant and product topologies, in terms of
atom properties changed, old interactions broken and new interactions formed. Updating existing interactions is
not possible, they have to be broken and re-created.

The directives that are valid in ``[moleculetype]`` are generally valid to use inside reaction templates. With the
notable difference in atom naming.


::

   [bonds]
   1:oh 2:oh 1 0.4 1000  ; a new harmonic bond formed in the reaction

   [angles]
   1:tail 1:oh 2:oh 1 85 100
   2:tail 2:oh 1:oh 1 85 100

Atom properties can be changed using the ``[redefine]`` directive. In each line of this directive,
a keyword specifies what is changed (name, type, mass, charge), followed by the new value.

::

   [redefine]
   1:reactive name BRD type SN1
   2:reactive name BRD type SN1

Existing interactions can be removed too, using the ``[break]`` and ``[update]`` directives.
Each line in these directives contains a list of atoms, in which the interactions should be broken
according to the rules of the directives.

::

   ; for a different reaction containing two atoms, left and right
   [break]
   1:left 1:right

The break directive completely disconnects the molecule between two
specified atoms. The bonds, angles, dihedrals, exclusions that contain
both atoms will all be removed.

The update directive has the same form as the break directive, but is more careful. It only breaks interactions
that only contains the specified atoms. So, a update directive with only two atoms specified will break the bond and
exclusion, but not angles or dihedrals that contain them and other atoms as well.

Lastly, reaction templates should be tested. It’s recommended to run reactive simulations with the
:doc:`/autoapi/martini_daemon/ReactionReporter`, as it reports each reaction, along with the frame it happens,
the internal fragment ID of reactants, and a list of atoms within the fragment.
For simple setups, the modification templates can be debugged using :doc:`/autoapi/martini_daemon/SystemDump`. It is
recommended to do so first in a small-scale simulation of a single reaction.

Reporters
---------

To produce output from simulations, that can be visualized or analyzed, reporters have to be added to simulations.
Reporters are Python class instances that inherit the ``Reporter`` base class, and their methods get called by
``Simulation`` at specified events during simulations.

Reporters can be broadly divided into distinct categories:

* Some perform reporting at a pre-defined trajectory frequency:
    * This frequency is the ``traj_frequency`` argument to the constructor of ``Simulation``.
    * This frequency is the same per-simulation to make analysis easier -- frame indices are synchronized across reporters.
    * For a successful simulation, there is a frame written at the start, when the current MD step modulo frequency is 0, and at the very end of the simulation.
    * :doc:`/autoapi/martini_daemon/VariablesReporter` reports thermodynamic variables for each trajectory frame.
    * :doc:`/autoapi/martini_daemon/TrajectoryReporter` reports atom positions and the periodic box for each trajectory frame. The format is implied from the file extension.
    * :doc:`/autoapi/martini_daemon/TopTrajReporter` creates a "topology trajectory", reporting atom properties and bonds for each trajectory frame. See :doc:`/toptraj` for details.
    * :doc:`/autoapi/martini_daemon/FragCountReporter` reports the number of fragments at each trajectory frame. It also adds the current total number of fragments to the interactive line on stdout.

* Some perform reporting related to reactions happening in the system:
    * :doc:`/autoapi/martini_daemon/LocalMinimizer` is not a traditional reporter. It locally minimizes the energy after reactions.
    * :doc:`/autoapi/martini_daemon/ReactionReporter` logs all reactions and reactants to a file as they happen.
    * :doc:`/autoapi/martini_daemon/ReactionEnergyReporter` reports thermodynamic variables before and after reactions. Optionally, it can write coordinates too, which can be helpful to debug local minimizations.

* Some of the reporters are there to help debug reaction templates:
    * :doc:`/autoapi/martini_daemon/SystemDump` prints all atom properties and a list of bonded forces in the system at the simulation start and after each frame with reactions. It is advised to use this on small systems only.
    * :doc:`/autoapi/martini_daemon/FragmentsDump` prints all fragments (successful graph matches).

Note, technically a reporter can do something at both trajectory frames and when reactions happen,
it is up to the implementation to choose which callbacks to hook on.

User-defined reporters are supported, they can use the public API of simulation
and all its public attributes to perform tasks during callbacks. Callbacks are called in the order
of the reporters list that is passed to Simulation. See :doc:`/extending` for more detail.

Analysis and Visualization
--------------------------

Standard tools can be used:

* For basic thermodynamic information the output of :doc:`/autoapi/martini_daemon/VariablesReporter` can be read as if
  it was a ``.csv`` file.
* ``Simulation.save_geometry()`` can be used to save a geometry snapshot (in a ``.gro`` or ``.xyz`` format).
* Standard ``.xtc`` or ``.trr`` files can be written using :doc:`/autoapi/martini_daemon/TrajectoryReporter`,
  which can be analyzed using GROMACS, mdtraj or MDAnalysis.
* If a ``.tpr`` file is needed, it is recommended to create a ``.tpr`` file using GROMACS, given the same topology,
  with any ``.mdp`` file. This should work alright for many analysis tools.

However, analyzing and visualizing based on the changes reactions introduce in the topology can be difficult with
pre-existing tools. See :doc:`/toptraj` for some of the ways Martini Daemon facilitates the analysis and visualization
of simulations with changing topologies.

For custom needs, a lot of information can be accessed through the Python API.
This can be accessed from the run script, or from callbacks in custom Reporters. See ``/extending``.

Checkpoint system
-----------------

TODO
