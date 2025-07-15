documentation:
- API documentation
- Reference (based on report submitted)
- User guide
- Some write up for the examples
- consistent style guide
- linting
- docstring everything + make pydoc work

testing overhaul:
- verify dihedral detection testing again
- Modification testing
- checkpoint testing

code quality issues:
- forces even more modular (return mm.Force, not set self._force_obj)
- reporter 2
- clean up helper
- adjust reaction rate counting to filter out when overlapping fragments react independently (currently rejected by the modification algorithm, right?)
- improve the passing around of values, remove all default values
	except in Simulation
- checkpoint reporter load() should become part of the initializer
- fix TODOs in code

T* breakup and rework:
- rearrange frag and rx info in T* to try to be always dense, benchmark change
- detection2 debug and benchmark
- make top parser fully internal, don't rely on it in tests
- break up T*, improve the QoL when using daemon only as a martini impl

Later:
- the full range of rate control functions
- don't reinit on only LJ type change / charge change (if possible)
- cleaner should-we-reinit logic
- relative rate control groups
- remove optional atoms if possible
- half reaction system
- multi reactant reactions cleaner with detection2
- friendlier angle/dihedral syntax
