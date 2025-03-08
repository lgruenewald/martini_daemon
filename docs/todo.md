- D/M unit testing
	- detection testing - single run of D algo
		- json recipe for list of reactions and list of atoms involved
	- modification testing - single run of D/M algo
		- fragment list similar to how graphs do it
- multiple reactions with the same starting materials, with equal standing (independent of their order)
- velocity rx_condition, better rx_condition syntax
- undoing reactions instead of [rx_break] and [rx_update]
- update documentation

Optimizations
- early quits in detection_over_types on more than 2 reactant situations when part of the conditions do not get fulfilled
- optimization: the same graphs per molecule type when doing it at the start
	- cache results for this and reuse
- do not add S* forces that are never going to change (reduced memory usage when dilute systems)
- S* should try to update forces instead of remove/adding when possible
  - allow this for vsites, constraints
- neighborlist for the D algorithm
- cythonize and parallelize bottlenecks

Correctness
- limiter should not depend on the order of atoms in the topology
- verification for the "equivalent" keyword
- reaction rates should be independent of how often the D/M algorithm is done

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
