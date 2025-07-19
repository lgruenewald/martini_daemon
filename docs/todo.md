documentation:
- Reference (based on report submitted)
- User guide
- Some write up for the examples

testing overhaul:
- verify dihedral detection testing again
- Modification testing
- checkpoint testing + integrate save()/load() into constructors

code quality issues:
- make a list of valid filters, reject unknown filters
- forces even more modular (return mm.Force, not set self._force_obj)
- reporter 2 - remove version 2 if not any better
- clean up helper
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
- the full range of rate control functions
- itp writer
- pbc whole
- don't reinit on only LJ type change / charge change (if possible)
- cleaner should-we-reinit logic
- relative rate control groups
- remove optional atoms if possible
- half reaction system
- multi reactant reactions cleaner with detection2
- friendlier angle/dihedral syntax
