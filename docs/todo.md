cythonization
1. index_pair solution for instantiate_over_existing
2. update topstar to use new rx templates and number only fragments
3. update detection.pyx
4. debug python <=> cython interface



helper folder:
all bond -> vmd bond
pbc whole
Benchmarking logs:
% slowdown of MD steps, and how long each MD step took (and their distribution)
D/M cost in milliseconds (most transparent and isolated way to do it)
average number of reactions per D/M step
% cost of each component in both reaction and non reaction case
% slowdown in total
temperature or other variables over time


- better quitting behavior in benchmarking
- linting, type checking
- everything should use "atoms", not "particles"

- optimization: the same graphs per molecule type when doing it at the start
	- cache results for this and reuse

- S* should try to update forces instead of remove/adding when possible
  - allow this for vsites, constraints

fixed reaction rates:
- remove limiter, probability, add energy barrier / velocity
- benchmark reaction rate vs D/M frequency, make sure it's constant at the default D/M frequency
- make sure energy barrier controls reaction rate completely
- reaction rates should be independent of how often the D/M algorithm is done within reasonable bounds

- D/M unit testing
	- detection testing - single run of D algo
		- json recipe for list of reactions and list of atoms involved
	- modification testing - single run of D/M algo
		- fragment list similar to how graphs do it
- multiple reactions with the same starting materials, with equal standing (independent of their order)
- velocity rx_condition
- update documentation

Correctness
- verification for the "equivalent" keyword

Simulation stability
- Langevin integrator friction and stability?
- how far can we build bonds/angles compared to equilibrium? rough guidelines?

Tooling
- write itp's at the end
- checkpoints for simulations
- .mdp parser

What would we need to get support for all .top files? Most to least important.
- fix vsite1
- position restraints
- dihedral type 4
