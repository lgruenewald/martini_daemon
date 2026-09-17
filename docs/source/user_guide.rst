User Guide
==========

Introduction to run scripts
---------------------------

Simulations in Martini Daemon are run using Python scripts. These contain calls to Martini Daemon’s API,
specifying simulation parameters and input files. A simple run script is provided
below, which can be adjusted to meet various needs.

::

   #!/usr/bin/env python3

   from martini_daemon import Simulation, VariablesReporter, TrajectoryReporter, ReactionReporter, TopTrajReporter
   import openmm as mm

   with Simulation(
       # path to GROMACS Topology
       topology="system.top",
       # path to Starting geometry
       geometry="system.gro",
       # output filenames will be prefixed with this name
       sim_name="out",
       # total MD steps
       md_steps=400000,
       # .xtc and other outputs write frequency, in MD steps
       traj_frequency=5000,
       # detection/modification algorithm frequency, in MD steps
       dm_frequency=250,
       # report thermodynamic variables, other reporters go here too...
       reporters=[
           # each reporter typically writes its own file with its own specific extension
           VariablesReporter(),
           TrajectoryReporter(),
           ReactionReporter(),
           TopTrajReporter(),
       ],
       # integrator with temperature coupling built-in
       integrator=mm.LangevinMiddleIntegrator(
           # temperature
           298 * mm.unit.kelvin,
           # collision frequency
           1.0 / mm.unit.picosecond,
           # integrator time step
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
   ) as sim:
       # minimize energy first, saving the minimized coordinates to min.gro
       sim.context.minimize_energy()
       sim.save_geometry("min.gro")
       # generate velocities at 298 K
       sim.context.generate_velocities(298)
       # run simulation for md_steps
       sim.simulate()
       # write final coordinates, velocities and box to a .gro file
       sim.save_geometry("out.gro")

See :doc:`/autoapi/martini_daemon/Simulation` for all available arguments, attributes, and methods.

Output files are written by Reporters. These reporters will own output handles. These are closed automatically
if using the ``with Simulation(...) as sim:`` syntax, or manually by calling ``Simulation.finish()``.
All output file paths are prefixed with the value passed to the argument ``sim_name``.
By default, no final or post-minimization geometry is written. The script should include ``sim.save_geometry(...)``
calls if this is desired.

Some notes regarding run.py files:

* D/M frequency of 250 is recommended by default. Lower frequencies can create
  a large performance overhead. Do note that the D/M frequency will impact
  reaction rates to an extent, depending on the nature of the system,
  the reaction and the conditions used.

* Langevin integrators are recommended at the moment, as they are the most tested.

* It is possible to do minimization and equilibration using Martini Daemon as well.
  A separate ``min+eq.py`` can be made for this purpose, and dm_frequency can be
  set to 0 to disable reactions.

Including reactions
-------------------

A typical workflow for adding a reaction consists of the following steps:

1. The desired reaction should be broken down to a mechanism of elementary steps.
2. For each mechanistic step, the Martini topologies of the reactant and product molecules should be obtained.
3. The difference between the two should be written down as a list of new interactions, as well as old interactions to break.
4. A graph for the reactant should be constructed. This lets Martini Daemon recognize which beads together form a
   reactive functional group (**fragment**).
5. A reaction template should be defined. This consists of reaction conditions and topology modifications.
6. Reaction templates should be fine-tuned to get the desired conditions, rates and products.
   Note, that the time scale of MD simulations may require compromises, such as unrealistically
   high rates. However, a detailed discussion of this is beyond the scope of this document.

This user guide focuses on steps 4 and 5. The end result of these steps is a reaction template (``.rx``) file, containing
graph definitions, reaction conditions and topology modifications.
These files should be included using ``#include`` into a ``.top`` file, similar to ``.itp`` files.

.top and .rx files
------------------

GROMACS ``.top`` and ``.itp`` files are documented in GROMACS documentation.
Martini Daemon's parser implements a large subset of them that is relevant for Martini simulations.

Some things to note about the preprocessor implementation:

* ``#define`` replacements are supported but currently cannot change the number of tokens, they always map
  one token to one token, even if they are empty. ``#define`` replacements can be chained.
* ``#ifdef`` and ``#ifndef`` are supported, but the closing ``#endif`` must be in the same file as the opening
  ``#ifdef``.
* ``""`` and ``<>`` are implemented identically for ``#include``, both search locally and in the include path currently.
* If no include path is specified, Martini Daemon will try to auto-detect a GROMACS installation, and include its
  force field directory.
* Each file can only be included using ``#include`` maximum once. An error is raised otherwise.
* Not supported: preprocessor macros, ``#if``.

A naming convention for Martini Daemon specific parts of topologies, is to put them in a file with the extension
``.rx``. These files should contain the ``[graph]`` and ``[reaction]`` directives for the system.
The file extension is arbitrary, as ``#include`` emulates the C preprocessor behavior.
``DAEMON`` is always defined for Martini Daemon, so ``.rx`` files should be included as:

::

    #ifdef DAEMON
    #include "path/to/file.rx"
    #endif

This ensures that the same ``.top`` file can be used with GROMACS and Martini Daemon without modification.

Fragments
---------

Before reactions can be defined, Martini Daemon must identify potential
reactants in the system. A **fragment** is a group of connected beads representing
a reactant, matching a pre-defined pattern. They can be thought of as a separate
*layer* on top of the simulation, containing only labels for various groups of atoms.
The complete list of fragments, updated throughout the simulation,
represents all potential reactants currently in the system.
This list is kept up to date using the graph matching algorithm.
Therefore, fragments are defined in input files as graphs, with node and edge
filters to match for.

The basic syntax for defining graphs is done using the ``[graph]`` directive.
Lines within this directive start with a keyword, followed by keyword-specific parameters.
All graphs must be given a name using the ``name`` keyword.
Beads to be matched (graph nodes) should be defined with the ``atom`` keyword,
followed by the graph node name, the name filter, and the type filter.
The graph node name is the name that will be used to refer to that bead in ``.rx`` input files.
The matched interactions (graph edges) can be generic (e.g. ``bond``),
or specific (e.g. ``harmonic_bond`` or ``morse_bond``).
All beads should be connected with such interactions to form a single connected graph.

Below is a simple example, where a bead representing an alcohol group (name: ROH, type: SP3) is matched.
Another bead is included for angle conditions and angle forces.
Its graph node name is tail, with a name match pattern that matches any bead name starting with C and any type.
A generic bond interaction is matched. Given this graph example, Martini Daemon would then find all beads with names
ROH and type SP3, that are bonded to any bead with a name starting with C and any type.
This bond can be any bond type listed in ``[bonds]``, or a constraint.

::

    [graph]
    name alc
    atom oh   ROH SP3
    atom tail C* *
    bond oh tail

The graph matching algorithm will keep an always up to date list of all graph matches for the detection algorithm.
This means that once the **contract** specified by the graph description above is broken, the corresponding
fragment is removed. This is done by matching all beads at the start of the simulation,
and then re-calculating matches for beads that participate in reactions, as well as their "neighbors", under
the assumption that only reactions change the topology of the system, and only topology changes can invalidate
graphs.

Interaction filters with more than two beads are also supported. Angles and dihedrals can be matched as well.
For example, ``angle oh bead2 bead3`` will match any angle consisting
of the three specified beads, in any order.

Graph matches are valid only if all specified beads and interactions are matched. If a single
bead or interaction is missing, the graph match is considered invalid. Every bead, as well as every interaction can
only be matched once.

Note, that there is currently no backtracking in the interaction matching (there is backtracking in node matching).
This only poses a problem when there are multiple interactions that can match the same interaction, such as
a simultaneous ``bond a b`` and ``harmonic_bond a b``, for a molecule that has a harmonic bond and a different
bond between ``a`` and ``b``. It is possible, that ``bond a b`` matches the harmonic bond, leading to no valid
graph matches because the harmonic bond filter is unfilled. A good rule of thumb is to connect the same group of beads
with distinct filters only.

An important property of the graph matching algorithm is that it tries to be exhaustive. It considers all
possible matches, and adds them all to the list of fragments. Two matches are considered identical if, for each
matched graph node, the same bead index within the system is matched. For our example, this means, that if there
are two different beads named C (or C1, CA, ...) bonded to the alcohol bead, two different matches can be obtained.
Both matches, by default, will be added to the fragment list as viable reactants that fulfill the contract specified
in the graph description. If this is undesirable, it is important to make graphs that do not have such ambiguity.
There are multiple strategies for this:

* Use stricter name and type filters than in this example. If possible, change the names of the beads in the .itp
  if they represent fundamentally different beads.
* Include other beads in the graph that can be used to disambiguate.
* Add connections (bond type 5) to the molecule, and match them using the ``connection`` filter.
  This bond type does not add any potential to it, it is used only as a unit of topological information.
  Note, that visualizations might draw a bond between the two particles, which may be undesirable.
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
On exchanging the two bead indices, the algorithm will not consider it a new match. Without the equivalent keyword,
two matches would be found, with tail and tail2 exchanged, and both would be separately added to the fragment list.
Auto-detection of equivalence is currently not possible, because when the graphs are defined,
it is unknown how reaction templates use them, so this must be specified and verified manually.

Removing reactant graphs on reactions
-------------------------------------

It is important to define reactant graphs and reactions in a manner, such that during a reaction, they no longer will
be valid reactants that can still react in the same reaction (unless this is explicitly desired). This is done
implicitly in Martini Daemon by specifying a graph that no longer matches after the reaction.

Here is a list of general strategies for achieving this:

* In reactions that remove bonds, this can be automatically achieved by specifying a required bond between the two
  bonded beads.
* In reactions that change bead names or types, this can also automatically be achieved with bead name and type filters,
  if they no longer match after the reaction.
* In other cases, a powerful tool can be the inclusion of forbidden beads in graph specifications. The rest of this
  subsection will give an introduction to them.

Forbidden beads are specified similarly to regular beads, but they start with the `atom!` keyword. These beads
function identically to regular beads, with one difference -- how to find out if a potential match is valid.
Regular beads must be matched. Forbidden beads must not be matched.
Consequently, interaction filters including forbidden beads are only a recipe on how to look for them.
Likewise, if a forbidden bead cannot be matched, the graph match remains valid.
Forbidden beads are enabled by the eagerness of the
graph matching algorithm. The graph matching algorithm only considers whether a partial match is valid or not, once no
more beads can be added to it based on the filters.

Note that there is no way to group forbidden beads together -- that is, rules such as "it is forbidden to have a
bead of this name and type that is bonded to this other bead with this name and type" are not possible. Another
way to think about this restriction is that forbidden beads can only represent a "one deep" layer around the normal
beads. This is important, so a graph can only become forbidden if a change occurred to its direct
neighbors. Therefore, graphs only need to be recalculated if their direct neighbors change. In order to
prevent users from attempting to group forbidden beads, there is an error message raised if there is an interaction
filter connecting a forbidden bead to another forbidden bead.

Note that forbidden beads can show up as a ``-1`` in some reporter outputs, since they are represented with a -1
internally in the fragment list.

Optional beads
--------------

Optional beads allow for some extra complexity in some cases. They function similarly to forbidden beads,
with the difference that the graph is valid regardless of whether they are there or not.

Rules for optional beads:

* Optional beads are defined with the ``atom?`` keyword, followed by graph node name, name filter, type filter.
* Missing optional beads appear as ``-1`` in some reporter outputs, as they are represented with a -1 internally.
* If a reaction condition references a missing optional bead, the condition is ignored.
* If a modification template entry references a missing optional bead, the whole entry is ignored.
  An exception to this is ``[update]``, which will proceed with the rest of the specified beads
  (this makes more sense for that specifically).
* Due to the eagerness of the graph matching algorithm, if an optional bead can be matched, only the graph match
  including the optional bead is added to the fragment list.
* Optional beads cannot be "grouped", that is interaction filters within the graph can only contain at most one
  optional bead per interaction filter. They also cannot be grouped with forbidden beads. This facilitates
  graph recalculation on direct neighbor change only.
* Due to neighbor recalculation, if an existing fragment with a missing optional bead suddenly has an optional bead
  available to it, it will be recalculated to include it.

Debugging graphs
----------------

The reporter :doc:`/autoapi/martini_daemon/FragmentReporter` has a detailed option, designed to facilitate debugging
the graph matching algorithm, by printing the list of matched fragments.

Graph syntax reference
----------------------

All lines in this directive should start with a keyword, followed by parameters. Below is a list of valid keywords.

- ``name graph_name`` - specify the graph name

- ``atom graph_node_name name_filter type_filter`` - mandatory beads
- ``atom? graph_node_name name_filter type_filter`` - optional beads
- ``atom! graph_node_name name_filter type_filter`` - forbidden beads

Graph node names can contain ASCII letters, numbers and underscores.

Name and type filters support the following characters or wildcards:

* alphanumeric characters and underscores will match themselves only
* ``*`` for any strings (including empty strings)
* ``?`` for any single character
* ``{123}`` will match one of the characters included inside the braces.

- ``equivalent graph_node_name1 graph_node_name2 ...`` - specify a group of beads, which when
  exchanged, do not represent a different graph match

- ``<interaction_name> graph_node_name1 graph_node_name2 ...`` - specify a group of beads,
  which are connected by the interaction filter. See section about
  interaction filters for valid filters.

The valid interaction filters are given below. Note, that User-defined forces can also be matched. This list is for
the built-in forces only.

- ``bond``: any bond listed here
- ``harmonic_bond``: bond type 1 and 6
- ``g96_bond``: bond type 2
- ``morse_bond``: bond type 3
- ``cubic_bond``: bond type 4
- ``connection``: bond type 5
- ``fene_bond``: bond type 7
- ``distance_restraint``: bond type 10
- ``constraint``: constraint type 1 and 2
- ``angle``: any angle listed here
- ``harmonic_angle``: angle type 1
- ``g96_angle``: angle type 2
- ``cross_bond_bond``: angle type 3
- ``cross_bond_angle``: angle type 4
- ``urey_bradley``: angle type 5
- ``quartic_angle``: angle type 6
- ``linear_angle``: angle type 9
- ``restricted_angle``: angle type 10
- ``dihedral``: any dihedral listed here
- ``proper_dihedral``: dihedral type 1, 4 and 9
- ``improper_dihedral``: dihedral type 2
- ``rb_torsion``: dihedral type 3
- ``fourier_dihedral``: dihedral type 5
- ``restricted_dihedral``: dihedral type 10
- ``combined_bending_torsion``: dihedral type 11
- ``virtual_site``, ``vsite``: any virtual site, matches the virtual site as well as constructing beads
- ``vsite1``: virtual_sites1 type 1
- ``vsite2``: virtual_sites2 type 1
- ``2fd``: virtual_sites2 type 2
- ``vsite3``: virtual_sites3 type 1
- ``3fd``: virtual_sites3 type 2
- ``3fad``: virtual_sites3 type 3
- ``3out``: virtual_sites3 type 4
- ``4fdn``: virtual_sites4 type 2
- ``com``, ``center_of_mass``: virtual_sitesn type 2
- ``weighted_average``: virtual_sitesn type 1 and 3
- ``pair``: only matches pairs
- ``cmap``: only matches cmap
- ``exclusion``: only matches exclusions

Reaction templates
------------------

Reaction templates define:

* Reactants that should react,
* Reaction conditions for the reaction algorithm,
* Topology modifications to perform.

All reaction templates should start with the ``[reaction]`` directive.
This directive must contain a single line, specifying a unique reaction name.
Reaction names must not conflict with existing molecule type names, or other reaction names.
Reaction names should only contain alphanumeric characters and underscores.

::

   [reaction]
   example_reaction




The list of reactants should be given as a list of graph names in the ``[reactants]`` directive.
This directive must come first after each ``[reaction]`` directive. This directive is mandatory for each reaction.
It must contain a single line, with space separated graph names.
The fragments in the fragment list with this name will be considered during the detection algorithm.
Up to three reactants are supported. All reactant graphs must be defined earlier in the input file.

::

   [reactants]
   alc alc

Reaction conditions
-------------------

Reaction conditions are specified in the ``[conditions]`` directive.
This directive must come after the ``[reactants]`` directive. This directive is mandatory for each reaction.
Each line in this directive specifies a single condition, starting with a keyword, followed by parameters.

Graph nodes in these conditions are specified using the ``reactant_index:graph_node_name`` syntax. Reactant index
should be a 1-based index of the reactant whose bead is being referenced. Graph node name should be the name
as specified by the second token in ``atom`` keyword lines. The reactant index and node name must be
separated by a colon, with no whitespace between. Referencing optional atoms is possible.
If an optional atom is missing from the graph, that condition is ignored. Referencing forbidden atoms is not allowed,
as they are never present in any fragment.

An example set of conditions can be found below.

::

   [conditions]
   r_max 1:oh 2:oh 0.5
   r_min 1:oh 2:oh 0.3
   angle 1:tail 1:oh 2:oh 70 to 100
   angle 2:tail 2:oh 1:oh 70 to 100
   dihedral 1:tail 1:oh 2:oh 2:tail 50 to 100
   p 0.5

The possible reaction conditions are:

* ``r_max bead1 bead2 distance``

    * reactions above the specified maximum distance (in nanometers) will be rejected

* ``r_min bead1 bead2 distance``

    * reactions below the specified minimum distance (in nanometers) will be rejected

* ``angle bead1 bead2 bead3 min to max or min2 to max2``

    * Only angles between ``min`` and ``max`` are allowed.

    * ``min`` and ``max`` should be in degrees, and must be between 0 and 180 (inclusive).

    * Multiple allowed ranges can be specified using the ``or`` keyword between them.

    * If multiple ranges are specified, they cannot overlap.

    * For a set of beads, only one ``angle`` condition is allowed. All ranges for these beads must be specified on
      a single line.

* ``dihedral bead1 bead2 bead3 bead4 min to max (or min2 to max2)``

    * Identical syntax to angles, but for dihedrals.

    * Any values can be used. They will be put back in the -180 to 180 range by the parser.

    * If the larger value comes first, it will be assumed that it should "wrap" around the period, that is specifying
      170 to -170 will mean that anything "after" 170 and "before" -170 is allowed. This allows 170 to 180 and
      -180 to -170 in practice.

* ``p probability``

    * adds a random probability of accepting the reaction, which is evaluated after all other checks
      have been met. Should be between 0 and 1.

	* This probability condition is rolled independently at each attempt independently. Reaction rates will still
	  be affected by the Detection/Modification frequency and geometry conditions.

All geometry conditions take the periodic boundary condition into account,
calculating them based on minimum-distance images.

A maximum distance condition ``r_max`` is required to "connect" all reactants, if there are multiple.

Arbitrary reactant orientation conditions can be constructed from a combination of these rules.
The following rule should be kept in mind when designing them:
each condition is evaluated separately, and a single failing condition will reject the reaction.

Note that a reaction can also be rejected for other reasons. Overlapping (that is they share a bead) fragments
cannot react with eachother. During the modification algorithm, all overlapping fragments will get invalidated as part
of the modification algorithm, before running the graph matching again.
This invalidation means that overlapping fragments cannot participate in separate reactions in a single D/M step,
nor can a single fragment participate in two different reactions. This is intentional, so that each reaction guarantees
that the contract (graph) of its reactants was upheld and had no undesired changes introduced in the same D/M step.

Note that there are no guarantees made for the order in which reactions are checked, neither of determinism, nor of randomness.

The following aspects should be considered when making reaction conditions:

* If forming a bond, generally there should be a maximum distance requirement for the given pair of beads.

* Forming a bond far away from the equilibrium bond length may result in a too high potential energy right after
  the reaction, or may result in other molecules still being "sandwiched" between the two reactants.

* If removing a bond, it can be useful to have a minimum distance requirement for the given pair of beads.

* If forming an angle, sometimes it is necessary to have an angle condition.

* Sometimes, to get the right product geometry, distance, angle or dihedral reaction conditions are needed.

* Some multi-functional molecules may self-react in undesirable ways if such conditions are too broad.

* When multiple reactants are specified, a maximum distance
  must be present to "connect" them. This is because the detection algorithm constructs a KdTree of all fragment
  locations in space, which is used to accelerate it. Without a maximum distance condition, all combinations
  of reactants would need to be considered.

Note that new reaction conditions can be added to this package in the future. Please open an Issue on Github.

Reaction rates
--------------

When tuning the reaction rates, be aware that:

* Not only reaction rates relative to other reactions have an effect on the simulation results,
  but also relative to all processes happening in the simulation, including diffusion.

* Different geometry conditions will result in different reaction rates. It may be useful to set reaction conditions
  at geometries that are somewhat higher in the potential energy surface, emulating a sort of "activation energy".

* Strict geometry conditions reduce reaction rates.
  Generally angle conditions have a larger impact than distance conditions,
  while dihedral conditions have an even larger impact.

* If using the probability reaction condition, the same reaction probability may result in different rates based on
  the other conditions, as the probability is only applied "on top" of any other conditions. A lower probability
  reaction can still be faster, depending on the concentration and geometry conditions.

Modification templates
----------------------

Reactions should contain a list of changes needed to be made to the topology after such a reaction occurs.
Currently, this is a manual
process. Users should first write down the difference between the reactant and product topologies, in terms of
bead types changed, old interactions broken, and new interactions formed. Updating existing interactions is
not possible, they have to be broken and re-created.

The directives that are valid in ``[moleculetype]`` are generally valid to use in ``[reaction]`` directives.

The important differences are:

* Beads are indexed using the ``reactant_index:graph_node_name`` syntax, identically to reaction conditions.
* Interactions referencing missing optional beads will be skipped.
* Adding virtual sites and constraints is not supported.
* The ``[atoms]`` directive is not supported. See ``[redefine]`` below on how to change bead properties.

::

   [bonds]
   1:oh 2:oh 1 0.4 1000  ; a new harmonic bond formed in the reaction

   [angles]
   1:tail 1:oh 2:oh 1 85 100
   2:tail 2:oh 1:oh 1 85 100

A note on exclusions.
Bonds added during reactions will also generate exclusions, just like in ``[moleculetype]``.
Exclusions can also be added manually using the ``[exclusions]``.
Currently ``nrexcl`` other than 1 is not supported in Martini Daemon.

Redefine
--------

Bead properties can be changed using the ``[redefine]`` directive.
Each line in the ``[redefine]`` directive starts with specifying which bead to change, specified using the
``reactant_index:graph_node_name`` syntax.
The rest of the line contains key value pairs, separated by whitespace.
The keys ``name``, ``type``, ``charge`` and ``mass`` are supported.

Changing the type alone will not affect the charge or mass! Even if the default charge/mass of the original
type was used. Mass and charge need to be explicitly changed, if changing them is desired.

::

   [redefine]
   1:reactive name BRD type SN1
   2:reactive name BRD type SN1

Break and Update
----------------

Existing interactions can be removed too, using the ``[break]`` and ``[update]`` directives.
Each line in these directives contains a list of beads, in which the interactions should be broken
according to the rules of the directives. Each of these lines describes a **break group** or **update group**.

::

   ; for a different reaction, let's say the reactant contains two bonded beads, left and right
   [break]
   1:left 1:right

During the modification algorithm, each break and update group is executed one by one. They choose interactions
that include those beads and remove them. The difference lies in the rules which govern the choice of these
interactions:

* **break** groups remove interactions which contain **all** beads in the group,
* **update** groups remove interactions which contain **only** beads from the group.

The naming is based on their intended use case:

* break is intended to be used to completely disconnect a set of particles (usually a set of 2). All bonds, angles and
  dihedrals will be broken. Note, that it doesn't traverse the whole interaction graph, so the two sides of
  the molecule can still remain connected through other beads not included in the break group.
* update is intended to be used to change interactions between a set of beads. If two beads are specified, it will
  only break the various bonds and exclusions between them, not the angles or dihedrals that also include other
  beads. This way, new bonds can be added to replace the old ones.

Optional beads are supported. A missing optional bead will:

* Disable the break group, as no interaction will contain beads that are not there.
* Still allow the remaining beads to be processed as an update group.


Debugging reaction templates
----------------------------

Reaction templates should be tested. It is recommended to run reactive simulations with the
:doc:`/autoapi/martini_daemon/ReactionReporter`, as it reports each reaction, along with the frame it happens,
the internal fragment ID of reactants, and a list of atoms within the fragment.
For simple setups, modification templates can be debugged using :doc:`/autoapi/martini_daemon/SystemDump`. It is
recommended to do so first in a small-scale simulation of a single reaction.

Reporters
---------

To produce output from simulations, that can be visualized or analyzed, reporters have to be added to simulations.
Reporters are Python class instances that inherit the ``Reporter`` base class, and their methods get called by
``Simulation`` at specified events during simulations.

Reporters can be broadly divided into distinct categories:

* Some perform reporting at a predefined trajectory frequency:
    * This frequency is the ``traj_frequency`` argument to the constructor of ``Simulation``.
    * This frequency is the same per-simulation to make analysis easier -- frame indices are synchronized across reporters.
    * For a successful simulation, frames are written at the start, when the current MD step modulo frequency is 0, and at the very end of the simulation.
    * :doc:`/autoapi/martini_daemon/VariablesReporter` reports thermodynamic variables for each trajectory frame.
    * :doc:`/autoapi/martini_daemon/TrajectoryReporter` reports atom positions and the periodic box for each trajectory frame. The format is implied from the file extension.
    * :doc:`/autoapi/martini_daemon/TopTrajReporter` creates a "topology trajectory", reporting atom properties and bonds for each trajectory frame. See :doc:`/toptraj` for details.
    * :doc:`/autoapi/martini_daemon/FragmentReporter` reports the number of fragments at each trajectory frame. It also adds the current total number of fragments to the interactive line on stdout.

* Some perform reporting related to reactions happening in the system:
    * :doc:`/autoapi/martini_daemon/LocalMinimizer` is not a traditional reporter. It locally (reacting atoms and their environment) minimizes the energy after reactions.
    * :doc:`/autoapi/martini_daemon/GlobalMinimizer` is also a minimizer, but it runs OpenMM's minimization algorithm on the entire system. This should be considered experimental.
    * :doc:`/autoapi/martini_daemon/GlobalIntegratorMinimizer` is a minimizer, which runs any OpenMM integrator for a few steps after each reaction. These extra steps are not counted toward simulation time. This should be considered experimental.
    * :doc:`/autoapi/martini_daemon/ReactionReporter` logs all reactions and reactants to a file as they occur. 
    * :doc:`/autoapi/martini_daemon/ReactionEnergyReporter` reports thermodynamic variables before and after reactions. Optionally, it can write coordinates too, which can be helpful to debug local minimizations.

* Some of the reporters are there to help debug reaction templates:
    * :doc:`/autoapi/martini_daemon/SystemDump` prints all atom properties and a list of bonded forces in the system at the simulation start and after each frame with reactions. It is advised to use this on small systems only.

Note, technically a reporter can do something at both trajectory frames and when reactions happen,
it is up to the implementation to choose which callbacks to hook onto.

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

For custom needs, more information can be accessed through the Python API.
This can be accessed from the run script, or from callbacks in custom Reporters. See ``/extending``.

Checkpoint system
-----------------

Martini Daemon includes a checkpoint system, based around :doc:`/autoapi/martini_daemon/CheckpointReporter`.
This reporter writes Martini Daemon specific ``.chk`` files, which contain (some of) the simulation state.
Checkpoints currently should be viewed as subject to change, as they currently only store particle positions and
velocities, the periodic box, and some simulation metadata. The topology is restored by loading the original ``.top``
file, and *replaying* all reactions that happened, based on the output of
:doc:`/autoapi/martini_daemon/ReactionReporter`.

When continuing simulations, writing will continue by appending to existing files. If there is extra output in the
output files since the checkpoint was made, the output files are truncated to provide a continuous output file.
The user should ensure proper backup practices are followed before continuing the simulation.

Loading checkpoints should happen using :doc:`/autoapi/martini_daemon/CheckpointLoader`. The first argument to the
constructor should be the path to the checkpoint file. All other arguments will be forwarded to the constructor
of :doc:`/autoapi/martini_daemon/Simulation` and should be identical to how the simulation ran before the
checkpoint was made.

::

   #!/usr/bin/env python3

   from martini_daemon import Simulation, VariablesReporter, TrajectoryReporter, ReactionReporter, TopTrajReporter, CheckpointLoader
   import openmm as mm

   with CheckpointLoader(
       "out.chk",
       # path to GROMACS Topology
       topology="system.top",
       # path to Starting geometry
       geometry="system.gro",
       ...


Energy minimization after reactions
-----------------------------------

Topology changes after reactions often create a sharp change in potential energy and forces in the system.
This is usually mitigated by the temperature coupling scheme, but in some cases this is not enough. In certain cases,
numerical instability can lead to crashes related to this. Martini Daemon features a :doc:`/autoapi/martini_daemon/LocalMinimizer`,
designed to perform a short energy minimization only on the reacting beads and their immediate surroundings.
There are a few options exposed, which can be used to tweak this process:

* The choice of minimizer and its properties. Currently only :doc:`/autoapi/martini_daemon/LocalGradientDescent` is available.

* The (maximum) number of steps to perform.

* The choice of which beads are *movable* during the minimization.

    * By default, reacting beads and their immediate neighbors are movable.

    * A radius (in nanometers) can be specified.

    * The bond graph can be recursively traversed to mark the entire molecule that reacts movable.

* Constraints can be temporarily replaced with stiff harmonic bonds.

Below is an example with all these options specified.

::

    LocalMinimizer(
        LocalGradientDescent(
            initial_step_size_nm=0.01,
            etol=0.001,
            smoothing_factor=0.1
        ),
        minimization_steps=500, # maximum, if delta E < etol it can stop early
        r_movable=1., # radius for movability, in nm
        whole_molecule=True, # traverse bond network to make all movable?
        harmonic_constraints=True, # change constraints to harmonic bonds temporarily?
    )

Additionally, the ``[soft_core]`` directive can also be used within reaction templates to temporarily
turn on soft-core potentials during the energy minimization phase.
Each line within this directive specifies the beads ``reactant_index:graph_node_name``, the soft core lambda
and alpha parameters.

It is recommended to use soft core if NaN exceptions happen when an exclusion is removed during a reaction,
and local minimization alone does not help. A lambda parameter between 0.5 (softest) and 1.0 (no soft core whatsoever)
can be used. An alpha parameter of 0.5 is recommended (though in principle, slightly larger alpha also makes
the soft core softer, assuming lambda is not 1.0). The soft core formula is inspired by the soft core interactions
in GROMACS used for free-energy interactions
(https://manual.gromacs.org/current/reference-manual/functions/free-energy-interactions.html#soft-core-interactions-beutler-et-al).

Example:

::

    [soft_core]
    2:1 0.85 0.5
    2:5 0.85 0.5
    1:1 0.85 0.5
