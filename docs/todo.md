- adjust reaction rate counting to filter out when overlapping fragments react independently (currently rejected by the modification algorithm, right?)
- more Detection testing
- Modification testing
- checkpoint testing
- detection2 debug and benchmark
- rearrange frag and rx info in T*
- make top parser fully internal, don't rely on it in tests
- improve the passing around of values, remove all default values
	except in Simulation
- make a nice API that can be documented with pydoc
- forces even more modular (return mm.Force, not set self._force_obj)
- don't reinit on only LJ type change / charge change (if possible)
- cleaner should-we-reinit logic
- checkpoint reporter load() should become part of the initializer
- friendlier angle/dihedral syntax
- remove optional atoms if possible
- fix failing pairs test
- dihedral type 4
- cmap
- tabulated bonds
- break up T*, improve the QoL when using daemon only as a martini impl
- fix TODOs in code
- consistent style guide (especially regarding type hints)
- linting
- docstring everything

new features?
- half reaction system
- multi reactant reactions using neighbor list union
