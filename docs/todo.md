- poly systems
	- hyperbranched (with rates) - kinetic analysis
	- polyurethane + loop analysis + kinetic analysis
	- step polymerization - polystyrene?
- benchmark reaction rate vs D/M frequency

- rate systems:
	- decay
	- acid/base

- silicate systems:
	- silica with chains
	- tetrahedral silica - larger system, run longer, some extra analysis?

- bonus system:
	- mietini

DOCS:
- new user's guide - basic concept explanation + init+propagation+termination poly system
- how-to-guide - focus on analysis! - polyurethane + loop analysis + kinetic analysis
- how-to-guide - focus on graph! - tetrahedral network
- reference - explanation, architecture, public API reference

- D/M unit testing
	- detection testing - single run of D algo
		- json recipe for list of reactions and list of atoms involved
	- modification testing - single run of D/M algo
		- fragment list similar to how graphs do it

- rates testing

- sim stability benchmark

- reactions vs dm freq

- api change - less work, more exposed / transparency, split up constructor
	- pick n threads, gpu

- load() -> optional constructor argument

- move stuff to core and clean up the package API

- bug - what if the same interaction gets matched multiple times in the graph
- alternative friendlier dihedral/angle syntax ?

- optimization: the same graphs per molecule type when doing it at the start
	- cache results for this and reuse

- detection2?

- don't reinitialize if only LJ type and/or charge change (unpaired, unexcluded)

- fix pairs
- dihedral type 4
- modularize nonbonded
- finish other helpers
