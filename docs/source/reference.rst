Reference Guide
===============

This document strives to exhaustively describe the main file formats implemented in Martini Daemon.

* Some familiarity with Martini Daemon is assumed for this document.
* API reference can be found at :doc:``/autoapi/martini_daemon``
* Algorithms, such as the graph matching algorithm are also described in the User guide at :doc:``/user_guide``.

.top and .rx files
------------------

Gromacs ``.top`` and ``.itp`` files are documented in GROMACS documentation.

Martini Daemon's parser implements a large subset of them that is relevant for Martini Simulations.

Some things to note about the preprocessor implementation:

* ``#define`` replacements are supported, but they currently cannot change the number of tokens, they always map
    one token to one token, even if they are empty. ``#define`` replacements can be chained.
* ``#ifdef`` and ``#ifndef`` are supported, but the closing ``#endif`` must be in the same file as the opening
    ``#ifdef``.
* ``""`` and ``<>`` is implemented identically for ``#include``, both search locally and in the include path currently.
* If no include path is specified, Martini Daemon will try to auto-detect a GROMACS installation, and include its
    force field directory.
* Each file can only be ``#include``ed maximum once. An error is raised otherwise.
* Not supported: preprocessor macros, ``#if``.

A naming convention for Martini-Daemon specific parts of topologies, is to put them in a file with the extension
``.rx``. These files should contain the ``[graph]`` and ``[reaction]`` directives for the system.
The file extension is arbitrary, as ``#include`` emulates the C preprocessor behavior.
The define ``DAEMON`` is always defined for Martini Daemon, so ``.rx`` files should be included as:

::

    #ifdef DAEMON
    #include "path/to/file.rx"
    #endif

This ensures, that the same ``.top`` file can be used with GROMACS and Martini Daemon without modification.
Alternatively, if reactions should be toggled on/off in separate simulation, a different define, such as REACT
can be used. Such defines can be passed in the defines parameter of :doc:`/autoapi/martini_daemon/Simulation`.


Graphs
------

The ``graph`` directive contains the following possible lines, each
starting with a specific keyword:

- ``name graph_name`` - specify the graph name
- ``atom graph_node_name name_filter type_filter`` - mandatory atoms
- ``atom? graph_node_name name_filter type_filter`` - optional atoms
- ``atom! graph_node_name name_filter type_filter`` - forbidden atoms
- ``equivalent name1 name2 ...`` - specify a group of atoms, which when
  exchanged, do not represent a different graph match
- ``<interaction_name> name1 name2 ...`` - specify a group of atoms,
  which are connected by the interaction filter. See section about
  interaction filters for valid filters.

Graph Node Names can contain ASCII letters, numbers and underscores.

Name and type filters can contain the following substrings:

* alphanumeric characters and underscores will match specific strings
* ``*`` for any strings (including empty strings)
* ``?`` for any single character
* ``{123}`` braces will match one of the characters included within them.

Interaction filters
-------------------

- ``bond`` - any bond listed here

- ``harmonic_bond`` - bond type 1 and 6

- ``g96_bond`` - bond type 2

- ``morse_bond`` - bond type 3

- ``cubic_bond`` - bond type 4

- ``connection`` - bond type 5

- ``fene_bond`` - bond type 7

- ``distance_restraint`` - bond type 10

- ``constraint`` - constraint type 1 and 2

- ``angle`` - any angle listed here

- ``harmonic_angle`` - angle type 1

- ``g96_angle`` - angle type 2

- ``cross_bond_bond`` - angle type 3

- ``cross_bond_angle`` - angle type 4

- ``urey_bradley`` - angle type 5

- ``quartic_angle`` - angle type 6

- ``linear_angle`` - angle type 9

- ``restricted_angle`` - angle type 10

- ``dihedral`` - any dihedral listed here

- ``proper_dihedral`` - dihedral type 1, 4 and 9

- ``improper_dihedral`` - dihedral type 2

- ``rb_torsion`` - dihedral type 3

- ``fourier_dihedral`` - dihedral type 5

- ``restricted_dihedral`` - dihedral type 10

- ``combined_bending_torsion`` - dihedral type 11

- ``virtual_site``, ``vsite`` - any virtual site, matches the virtual site as well as constructing atoms

- ``vsite1`` - virtual_sites1 type 1

- ``vsite2`` - virtual_sites2 type 1

- ``2fd`` - virtual_sites2 type 2

- ``vsite3`` - virtual_sites3 type 1

- ``3fd`` - virtual_sites3 type 2

- ``3fad`` - virtual_sites3 type 3

- ``3out`` - virtual_sites3 type 4

- ``4fdn`` - virtual_sites4 type 2

- ``com``, ``center_of_mass`` - virtual_sitesn type 2

- ``weighted_average`` - virtual_sitesn type 1 and 3

- ``pair`` - only matches pairs

- ``cmap`` - only matches cmap

- ``exclusion`` - only matches exclusions

Reaction templates
------------------

Reactions are defined using the ``[reaction]`` directive. This directive should only contain one line,
specifying a unique reaction name. Reaction names must not conflict with existing molecule type names. Reaction names
should only contain alphanumeric characters and underscores.

Reactants
---------

Reactants are defined using the ``[reactants]`` directive. This directive should come first after each ``[reaction]``
directive, and is mandatory. It should contain a single line, with space separated graph names. Currently, up to three
reactants can be specified. All reactant graphs must be defined before in the input file.

Reaction conditions
-------------------

The possible reaction conditions are:

* ``r_max atom1 atom2 distance`` - reactions above the maximum distance (in nm) will be rejected
* ``r_min atom1 atom2 distance`` - reactions below the minimum distance (in nm) will be rejected
* ``angle atom1 atom2 atom3 min to max or min2 to max2`` - only angles between min and max are allowed
    * min and max should be in degrees, and must be between 0 and 180 (inclusive).
    * multiple allowed ranges can be specified using the ``or`` keyword between them
    * if multiple ranges are specified, they cannot overlap
    * for a set of atoms, only one ``angle`` condition is allowed = all ranges must be specified on one line
- ``dihedral atom1 atom2 atom3 atom4 min to max (or min2 to max2)`` - identical syntax to angles, but for dihedrals
    * but! the allowed range is anything! it will be put back in the -180 to 180 range by the parser
    * if the larger value comes first, it will be assumed that it should "wrap" around the period, that is specifying
      170 to -170 will mean that anything "after" 170 and "before" -170 is allowed. This allows 170 to 180 and
      -180 to -170 in practice.
- ``p probability`` - adds a random probability of accepting the reaction, which is evaluated after all other checks
    have been met. Should be between 0 and 1.

Atoms in these conditions are specified using the ``reactant_index:graph_atom_name`` syntax. Reactant index
should be a 1-based index of the reactant whose atom is being referenced. Graph atom name should be the name
of the atom in the graph, as specified by the second token in lines with the ``atom`` keyword. The two must be separated
by a colon, with no whitespace between. Referencing optional atoms is possible.
If an optional atom is missing from the graph, that condition is ignored. Referencing forbidden atoms is not allowed,
as they are never present in any fragment (since they are only used to create conditions for graph matching).

A maximum distance condition is required to "connect" all reactants, if there are multiple.

Each condition is evaluated separately. A single failing condition will reject the reaction.

Note, that a reaction can also be rejected for other reasons. Overlapping (that is they share an atom) fragments
cannot react. During the modification algorithm, if one reaction invalidates the fragment of another reaction,
the other reaction will also be rejected.

Bonds, angles, dihedrals
------------------------

Directives that add bonds, angles, dihedrals and similar interactions within ``[moleculetype]``
are also supported within reaction conditions.

The important differences are:

* Atoms are indexed using the ``reactant_index:graph_atom_name`` syntax, identically to reaction conditions.
* Interactions referencing missing optional atoms will be skipped.
* Adding virtual sites and constraints is not supported.
* The ``[atoms]`` directive is not supported. See ``[redefine]`` below on how to change atom properties.

Breaking bonds
--------------

Two directives are available for breaking interactions, the ``[break]`` and ``[update]``. Both share syntax,
they both contain lines with a list of whitespace separated atoms, picked using the ``reactant_index:graph_atom_name``
syntax. Each line specifies a new break group or update group respectively.

During the modification algorithm, each break and update group is executed one by one. They choose interactions
that include those atoms and remove them. The difference lies in the rules which govern the choice of these
interactions:

* **break** groups remove interactions which contain **all** atoms in the group,
* **update** groups remove interactions which contain **only** atoms from the group.

The naming is based on their intended use case:

* break is intended to be used to completely disconnect a set of particles (usually a set of 2). All bonds, angles,
  dihedrals will be broken. Note, that it doesn't traverse the whole interaction graph, so the two sides of
  the molecule can still remain connected through other atoms not included in the break group.
* update is intended to be used to change interactions between a set of atoms. If two atoms are specified, it will
  only break the various bonds and exclusions between them, not the angles or dihedrals that also include other
  atoms. This way, new bonds can be added to replace the old ones.

Optional atoms are supported. A missing optional atom will:

* Disable the break group, as no interaction will contain atoms that are not there.
* Still allow the remaining atoms to be processed as an update group.

Changing atoms
--------------

The ``[redefine]`` directive is used to change atom properties (currently name, type, charge, mass).

Each line in the redefine directive starts with specifying which atom to change, specified using the
``reactant_index:graph_atom_name`` syntax. After this, pairs of keywords followed by the new value are expected.
The keywords ``name``, ``type``, ``charge`` and ``mass`` are supported.

Example:

::

    [redefine]
    1:react name R type TC5 charge 0

Soft-core potentials for minimization
-------------------------------------

An additional feature within Martini Daemon is energy minimization after a reaction. Sometimes, changing the non-bonded
parameters can help with this. The ``[soft_core]`` directive can be used within reaction templates to temporarily
turn on soft-core potentials during the energy minimization phase.

Each line within this directive specifies the atom ``reactant_index:graph_atom_name``, the soft core lambda
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