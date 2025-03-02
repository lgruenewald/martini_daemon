A sunday's worth of work

- test and write examples

- optimization: the same graphs per molecule type when doing it at the start
	- cache results for this and reuse

Fragment atoms to ordered maps

- atoms inside a fragment shouldn't be just ordered, it should be a hashmap
	with optional names
- indexing in reaction should be 2D and allow these names
- instantiate products over a chosen ordered list of atoms instead of only indexing

Adding bonds without making new molecule types

- 0 products should be valid
- allow defining new bonds, angles, etc. in reactions

More flexible templates
- updating particle names during reactions
- early quits in detection_over_types on more than 2 reactant situations when part of the conditions do not get fulfilled
- reaction name should be usable for [frag_from], frags should be spawnable as products
- [atoms] should be optional, allow specifying just the number of atoms
- simplified [atoms] with less columns
- multiple reactions with the same starting materials, with equal standing (independent of their order)
- velocity rx_condition, better rx_condition syntax
- rename rx_..., frag_... to more sane names
- undoing reactions instead of [rx_break] and [rx_update]
- update documentation
- build the sticky martini system with the more flexible templates

Optimizations
- do not add S* forces that are never going to change (reduced memory usage when dilute systems)
- S* should try to update forces instead of remove/adding when possible
  - allow this for vsites, constraints

Rate correctness
- limiter should not depend on the order of atoms in the topology
- reaction rates should be independent of how often the D/M algorithm is done

Behavior
- better rx_break, rx_update, constraits to updating atom types (only when rx_update?)

Testing
- Unit test D/M algorithm

Simulation stability
- Langevin integrator friction and stability?
- how far can we build bonds/angles compared to equilibrium? rough guidelines?

Tooling
- write itp's at the end
- checkpoints for simulations
- .mdp parser

What would we need to get support for all .top files? Most to least important.
- fix pairs - clue: 1-2-3-4 bond network and 1-4 pair is bad, when removing 2-3 bond it somehow works?
- fix vsite1
- position restraints
- dihedral type 4
- tabulated potentials
- pbc that are not 90 degree angle
- cmap, just because martini_openmm can do it
