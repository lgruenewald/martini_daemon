








Relationships:

- DaemonParser builds T* from a static Topology (as defined in .top)
- T* builds and updates S*
- S* builds and updates the MD backend
- the D/M algorithm is currently a distinct part of T*
- the Simulation helper provides the API for the internals for common use cases
  similar to the openmm.app module

## DaemonParser

A parser that builds T* and S* from files on the disk. It can parse .itp and
.top files designed for the MARTINI force field in GROMACS, with some
extensions to allow for defining the templates.

## S*

Short for System Star, the layer that connects T* to a MD backend.

List of particles, bonds, angles, ...

Main task: abstract away access to openmm internals. Only use the openmm API
from this file, never from any other class. In theory this means that it's
easier to port martini_daemon to either different engines than openmm, or to
replace which parts of the API are used. It should try to also optimize
the internals as much as possible (reinitialize sparingly, only when needed).

Later should be optimized to not duplicate in-memory information about
non reactive particles.

## T*

Short for Topology Star, the layer that connects a static Topology 
(as defined in the input files) with S*, and can change during runtime.

A list of fragment types and reaction templates.

Fragments can originate from:
- molecules (the [moleculetype] directive)
  - they can be *instantiated*, either at the start during the [molecules]
  directive, or during reactions as defined in [rx_products], during this
  the particles and interactions defined in it get added to the system
  (or if it is during a reaction, existing atoms are updated to match
  the particles, rather than added)
- explicit definitions of fragments (the [frag] directive), these are defined
  as a (possibly reordered) ordered list of particles in another *parent*
  fragment, and get instantiated when the parent fragment is instantiated,
  in a recursive manner.
  This construct/layer of abstraction gives flexibility to templates.
  These can be used to, for example:
    - to abstract a "functional group" from different molecules
    to be able to perform the same reaction,
    - to handle symmetry (by adding the same fragment in multiple different
    order of atoms to the same parent fragment),
    - to change the order of atoms when considered in the reaction, so the
    order aligns with how the product [moleculetype] is defined

T* always keeps track of a list of fragments that are capable of reacting.
Fragments that have no reactions defined do not get added to this list
after instantiation. T* also has a defrag_list, which contains the list
of fragments each atom is in.

Reaction templates contain the following information, all used by the D/M
algorithm:
- reaction name (the first line in the [rx] directive)
- reactants ([rx_reactants] directive), can be 1 or 2
- products ([rx_products] directive)
- conditions ([rx_conditions] directive)
- options to break interactions ([rx_break], [rx_update] directives)

## D/M algorithm

Short for Detection/Modification algorithm.

### Detection algorithm

Checks whether the conditions for a reaction are fulfilled. Below is a list
of conditions for a reaction to be possible. A reaction is only possible
if all conditions are fulfilled. Any condition being failed means the detection
algorithm rejects the reaction.

Fragments that overlap (share an atom) cannot react. This condition is present
for all reactions. This fact can be used practically, for example to prevent
self reaction of a two-functional monomer.

If a probability condition is defined, every time the detection algorithm
is ran it rolls a random number from 0.0 to 1.0, and if the probability cutoff
is smaller than this random number, the reaction is rejected.

If a limiter is defined, for every reaction step in the simulation, the
specific reaction type can only happen up to the defined amount. This can
be useful to possibly speed up the detection step. Currently setting this
might bias the system to prefer reactions on the molecules that come first
in topology.

At least one distance cutoff for every reaction is mandatory. There is both
r_max possible for maximum distance and r_min possible for minimum distance.

An angle blacklist can be defined, which rejects reactions if the angle
defined by any 3 atoms in the reactants falls between the defined minimum and
maximum.

A dihedral blacklist can be defined, which rejects reactions if the dihedral
defined by 4 atoms in the reactants falls between the defined minimum and
maximum cutoffs.

In the future more conditions can be introduced, for example the kinetic
energy of the particle is planned.

### Modification algorithm

The modification algorithm will first break interactions defined in
[rx_break] and [rx_update]. Interactions will break if:
- If all atoms on an [rx_break] line are part of an interaction, the
  interaction is removed. This models for example bond breakage.
  Notably, it will remove angles, dihedrals, etc. that might have other
  atoms included too. Adding only a single atom on an [rx_break] line
  will remove all interactions for the given atom.

- If all atoms in an interaction are present in a [rx_update] line,
  the interaction is removed. This notably only includes interactions
  within the specified list of atoms on one line, if the interaction contains
  other atoms, it will not be removed. Adding only a single atom on a
  [rx_update] line will only remove interactions that only contain that one
  atom (since single atom interactions do not exist currently, this will
  do nothing).

Then, the modification algorithm will update the types, charge, mass of all
atoms according to the product [moleculetype]. If the type is * in the 
[moleculetype], the type is not updated.

The product particles are 1-1 mapped to the reacting fragments particles,
in the same order (use a [frag] to reorder if needed). When a product is
instantiated, new [frag]'s will be also recursively instantiated, so
reaction products can still be reactive.

## Simulation user experience // Tooling

- reporters in martini_daemon.reporters can do various cool things:
  - the bond_reporter can report the list of bonds (currently only those that
  do not cross the pbc), separately in every frame.



This document describes the additional directives that can be used for
creating reaction templates within martini_daemon.

# General recommendations

Martini_daemon specific things, such as [frag], [rx] directives, and the
product [moleculetype] directives (and all the directives within these
"parent" directives) should go to a separate file, currently recommended
with the ".rx" extension, which should be included in the main topology file,
after the .itp files.

If the topology files are being parsed by DaemonTopFile, the define DAEMON is
set, so putting the .rx files within an #ifdef DAEMON is possible, if the same
topology file is used in other software.


# Defining fragments - the [frag] directive

The first line in the [frag] directive should define the fragment name. The
fragment name should be unique (there should be no other [frag] or 
[moleculetype] with the same name). Fragment names can only contain ASCII 
letters, underscores and numbers. The [frag_*] subdirectives should
follow the [frag] directive directly for clarity.

## The [frag_atoms] directive

Optional. Must be after a [frag] directive was already defined. It applies
to the last [frag] directive. It is only valid to define this once per
[frag] directive.

Defines a list of atom types and name patterns to match against when
instantiating this fragment. If there is a mismatch, an error is raised.
Can be used for debugging templates, or to limit their misuse.

## The [frag_from] directive

Must be after a [frag] directive was already defined. If [frag_atoms] is present
for that [frag] directive, must come after [frag_atoms]. It applies to the
last [frag] directive.

Every line in this directive defines how this fragment should be instantiated.
The first column is the name of the parent fragment, and subsequent columns
contain the particle id's in the parent fragment that get mapped to this
fragment's particles. The number of particles in a fragment should be the same
for all lines, and same with the number of particles in [frag_atoms], if it
is defined.

All [frag_from] directives should come after the [frag] and [moleculetype]
directives they reference.

# Defining reactions - the [rx] directive

The first line in the [rx] directive should give the name of the reaction.
Reaction names can only contain ASCII letters, underscores and numbers.
The reaction name is only used for debugging at the moment, but in principle
they should be unique.

All [rx] directives should come after the [frag] and [moleculetype]
directives they reference, but before [system].

[rx_reactants], [rx_products] and [rx_conditions] are consdiered mandatory.
The directives within a [rx] directive should follow the [rx] directive
directly for clarity. Any of the [rx_*] subdirectives apply to the last
defined [rx] directive. It is only valid to define the [rx_*] subdirectives
once per [rx] directive, and they can come in any order.

## The [rx_reactants] directive

Contains the list of reactants for this reaction on a single line.
Can be [moleculetype] or [frag]. The numbering of the atoms for the product
and conditions will be the same as in the reactants. There can be 1 or 2 
reactants, other numbers of reactants are not supported.

## The [rx_products] directive

Contains the list of products for this reaction. Currently only 1 product
is supported, and it is recommended to make a new [moleculetype] entry
tailored to this specific [rx]. More flexibility might come later.

## The [rx_conditions] directive

List of reaction conditions to fulfill. See `docs/architecture.md`, section
Detection algorithm for details. One condition per line. Particles are indexed
in the following manner: first they index in reactant 1's atoms, then, if
there is a second reactant further indices index in reactant 2's atoms. So,
for example, if both reactants contain 3 particles, indices 1, 2 and 3 index
the first reactant and indices 4, 5 and 6 index the second reactant.

The following conditions are supported:

- `r_max`, arguments `i`, `j`, `r`, which rejects the reaction if the distance
  between particles `i` and `j` is above `r`.
- `r_min`, arguments `i`, `j`, `r`, which rejects the reaction if the distance
  between particles `i` and `j` is below `r`.
- `banned_angle`, arguments `i`, `j`, `k`, `min`, `max`, which rejects the
  reaction if the angle `ijk` (central atom `j`) falls between `min` and `max`.
  Both `min` and `max` are in degrees and the valid results for the angle
  calculation range from 0 to 180 degrees.
- `banned_dihedral`, arguments `i`, `j`, `k`, `l`, `min`, `max`, which rejects
  the reaction if the dihedral `ijkl` falls in the range `min` to `max`. Both
  `min` and `max` are in degrees. Valid results for the dihedral calculation
  range from 0 to 360 degrees (360 not included). If `min` is larger than `max`,
  it is assumed that it was meant to wrap around (for example, if `min` is 270
  and `max` is 90, it disallows it 0 to 90 and 270 to 360).
- `p`, argument `probability`, which rejects the reaction based on random
  chance. The `probability` is between 0.0 and 1.0, and represents the chance to
  accept the reaction (1.0 = always accept, 0.0 = always reject).
- `limit`, argument `count`, which rejects each reaction after `count` of this
  reaction have been accepted in this current reaction step.

## The [rx_break] and [rx_update] directives

Every line within these directives represents a group of particle indices to
be included in a break group. Each such group is used during a modification
algorithm to remove some interactions from reactants. Indices in these
directives work identically to how they work in [rx_conditions].

The difference between the two, is how interactions that are broken are
selected, see the section in `docs/architecture.md` for this for more
detail.

The main use case for [rx_break] is to specify a pair of particles, and to
break the molecule apart into two between those two particles. Note: if
additional exclusions were added when two monomers were connected, 
besides the main new bonded interactions, remember to break all these
exclusions as well, not just the bonded interactions.

The main use case for [rx_update] is to change how a part of the molecule
is constructed. It will break interactions within the group, but will not
break interactions that contain particles outside the group.

# Defining reaction products

Currently, it is recommended to write new [moleculetype] entries that define
a reaction product in the .rx file, specifically tailored with the reaction
template in mind. This is, because it is convenient to only define the
interactions that get formed during the reaction, while keeping most of
the interactions in the reactants that do not need updating. If an interaction
needs updating, they can be broken in the reaction template and readded by
the product.

## [atoms] directive in reaction products

An important difference between defining reaction products and just regular
[moleculetype] molecules is, that the [atoms] directive in reaction products 
does not represent new particles to be added to the system. The atom type,
charge and mass given here change the type, charge and mass of existing
particles. Other things, such as atom name, residue name, are currently
unchanged.

There are additional reaction product specific changes to the [atoms] directive
syntax in daemon:
- `*` is a valid atom type, which signifies that the atom type in the reactant
  should remain unchanged during the reaction. Note, that atom types are not
  patterns, only the `*` is recognized as a special value.

# Using non-daemon [moleculetype]s as reaction products

It is possible to use a predefined martini molecule as a reaction product.

To do this:
- add all reactant atoms to a single line within the [rx_update] directive
  for all reactants, so that all interactions that existed are broken first.
- add the name of the molecule to the [rx_product] directive

Defining new [moleculetype] directives is generally preferred, though, and
might be necessary, if the reactants and product contain constraints or
virtual sites, which cannot be broken or formed during reactions (also 
discussed in `docs/architecture.md`).

# Pattern matching

Pattern matching within the context of daemon is defined here. Since it is
mentioned in different later sections, it is described in a separate section.

A pattern matches a word if they are identical, but a few special characters
in patterns are also possible: `*` matches any string, `?` matches any single
character and curly braces (`{ab}`) match any of the characters enclosed
within the curly braces. A special character enclosed in curly braces matches
that special character literally.

Example use cases:
- `*` will match anything
- `TC*` will match any string starting with TC
- `C?` will match any string that contains C and one other character
- `C{123}` will match `C1`, `C2` and `C3`

Note: currently the unix fnmatch function is used to power this, which is also
used to match filenames, so it should be familiar.
