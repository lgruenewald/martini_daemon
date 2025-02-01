# Daemon architecture

Elements:
- DaemonParser
- S* (SysStar)
- T* (TopStar)
- D/M algorithm
- Simulation helper

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
<!-- TODO -->

<!-- TODO: later allow reaction rates to be independent of how often reaction
steps are done. -->

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

All fragments that contain atoms that are in a [rx_update] block will
also be deleted, since [rx_update] is used to "restructure" a part of
the molecule, signifying that its functionality is changing.
<!-- (TODO: only update atom types, charge when rx_update-->

Then, the modification algorithm will update the types, charge, mass of all
atoms according to the product [moleculetype]. If the type is * in the 
[moleculetype], the type is not updated.

The product particles are 1-1 mapped to the reacting fragments particles,
in the same order (use a [frag] to reorder if needed). When a product is
instantiated, new [frag]'s will be also recursively instantiated, so
reaction products can still be reactive.

## Simulation user experience // Tooling

<!-- - TODO, write the final Topology to a file, along with the final coordinates
  to allow continuing simulations with reactions, additionally TODO
  simulations with checkpoints -->
- reporters in martini_daemon.reporters can do various cool things:
  - the bond_reporter can report the list of bonds (currently only those that
  do not cross the pbc), separately in every frame.

## Limitations

- the number of particles cannot be changed during reactions
- constraints and virtual sites cannot be created or removed during reactions
  - cannot be created because of periodic boundary conditions inside openmm
  - cannot be removed because it messes indexing up
  - currently the parameters for them also cannot be changed, this limitation
  could be lifted later
- fragments that overlap can't react, there is no plans to support this,
  as this can be useful for building templates
- no checking for duplicate exclusions created between two particles during
  reactions (yet?) - do not create exclusions between pairs of atoms that
  are already excluded, this limitation will probably be lifted later
- reaction constraints and product interactions specified by the user
  should lead to forces during bond formation that don't blow up the system,
  and it's the user's responsibility to ensure this
  (martini_daemon ultimately should do what the user tells it to do)
