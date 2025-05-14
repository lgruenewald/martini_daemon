DOCS:
- new user's guide - basic concept explanation, simple example
- how-to-guide
	- more advanced polymerization
	- more advanced silica based guide
	- with analysis!
- reference - explanation, architecture, public API reference

- warn user if no fragments
- bug - what if the same interaction gets matched multiple times in the graph
- alternative friendlier dihedral/angle syntax ?

- benchmark reaction rate vs D/M frequency, make sure it's constant with default settings

- graph bug - trying to match the same interaction multiple times

- T* reporters auto generate graphs that are interesting

mol_fragment -> molecule
graph_fragment -> graph
particles -> atoms

helper folder:
working with top trajectories
all bond -> vmd bond
pbc whole
finalize itp

- optimization: the same graphs per molecule type when doing it at the start
	- cache results for this and reuse

- S* should try to update forces instead of remove/adding when possible
  - allow this for vsites, constraints

- D/M unit testing
	- detection testing - single run of D algo
		- json recipe for list of reactions and list of atoms involved
	- modification testing - single run of D/M algo
		- fragment list similar to how graphs do it

What would we need to get support for all .top files? Most to least important.
- fix pairs
- position restraints
- dihedral type 4
