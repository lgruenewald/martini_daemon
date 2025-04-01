code quality:
- mypy type checker
- linter - pylint, flake8, black
- naming conventions
- atoms vs particles - stick to one name - atoms, rename all "part, particles"
- coverage.py for coverage
- scalene for benchmarking
- cython

potential bugs:
- bug in bond reporter?
- better quitting behavior in benchmarking

fixed reaction rates:
- remove limiter, probability, add energy barrier / velocity
- benchmark reaction rate vs D/M frequency, make sure it's constant at the default D/M frequency
- make sure energy barrier controls reaction rate completely

- D/M unit testing
	- detection testing - single run of D algo
		- json recipe for list of reactions and list of atoms involved
	- modification testing - single run of D/M algo
		- fragment list similar to how graphs do it
- multiple reactions with the same starting materials, with equal standing (independent of their order)
- velocity rx_condition, better rx_condition syntax
- update documentation

Optimizations
- early quits in detection_over_types on more than 2 reactant situations when part of the conditions do not get fulfilled
- optimization: the same graphs per molecule type when doing it at the start
	- cache results for this and reuse
- do not add S* forces that are never going to change (reduced memory usage when dilute systems)
- S* should try to update forces instead of remove/adding when possible
  - allow this for vsites, constraints
- cythonize and parallelize bottlenecks

Correctness
- limiter should not depend on the order of atoms in the topology
- verification for the "equivalent" keyword
- reaction rates should be independent of how often the D/M algorithm is done within reasonable bounds

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
