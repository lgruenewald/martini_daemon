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
- checkpoint load testing
- reporter testing
- helper testing

- api change - less work, more exposed / transparency, split up constructor
- pick n threads, gpu
- make everything use DaemonSimulation as API, rename DaemonSimulation to Daemon
- make TopParser internal
- make forces even more modular / non inheritance

- don't reinitialize if only LJ type and/or charge change (unpaired, unexcluded)

- checkpoints: load() -> optional constructor argument

- move stuff to core and clean up the package API

- alternative friendlier dihedral/angle syntax ?

- optimization: the same graphs per molecule type when doing it at the start
	- cache results for this and reuse

- detection2?

- fix pairs
- dihedral type 4
- cmap
- finish other helpers
