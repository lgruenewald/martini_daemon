documentation:
- Reference (based on report submitted)
- User guide
- Some write up for the examples

testing overhaul:
- test custom forces
- verify dihedral detection testing again
- Modification testing
- checkpoint testing + integrate save()/load() into constructors

code quality issues:
- move custom stuff to plugins / separate repos
- finish/clean up helpers
- pbc whole constraints + vsites? when loading .gro files, warn if not possible
- forces even more modular (return mm.Force, not set self._force_obj)
- adjust reaction rate counting to filter out when overlapping fragments react independently (currently rejected by the modification algorithm, right?)
- improve the passing around of values, remove all default values
	except in Simulation
- fix TODOs in code or convert them to todo's here

Later:
- fix cmap, pairs
- rearrange frag and rx info in T* to try to be always dense, benchmark change
- detection2 style detection
- break up T*, improve the QoL when using daemon only as a martini impl
- make top parser fully internal, don't rely on it in tests
- the full range of rate control functions, relative rate ccontrol groups
- itp writer
- pbc whole
- don't reinit on only LJ type change / charge change (if possible)
- cleaner should-we-reinit logic
- remove optional atoms if possible
- half reaction system
- multi reactant reactions cleaner with detection2
- friendlier angle/dihedral syntax
