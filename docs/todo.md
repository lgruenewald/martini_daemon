More flexible templates

- reoptimize D algo for n reactants

- allow constructing fragments based on graphs, even without explicitly mentioning
the parent fragment name
  - allow "negative" constraints in these graphs (e.g. not bonded to X, bonded to up to n instances of Y)
- indexing in products, conditions should be 2D
- fragments with variable number of atoms, out of index reads in reaction templates
  should be ignored, if they are marked with a ?


- interactions should be per atom, not in Fragments
- build the sticky martini system with the more flexible templates
- fix rx_update removal of overlapping fragments, add removal of overlapping fragments to rx_break
- reaction name should be usable for [frag_from], frags should be spawnable as products
- 0 products should be valid
- overlapping fragments reacting, but all overlaps need to be explicit - [rx_overlap]
- [atoms] should be optional, allow specifying just the number of atoms
- allow defining new bonds, angles, etc. in reactions
- simplified [atoms] with less columns
- instantiate products over a chosen ordered list of atoms instead of indexing
- multiple reactions with the same starting materials, with equal standing (independent of their order)
- velocity rx_condition, better rx_condition syntax

.rx validation
- try to match name and type match at "compile time" to catch errors early
- construct list of all LJ types at the start
- compile reactions into a set of minimum operations - don't add/remove if change is possible
- allow changing of vsite and constraint parameters this way
- error on attempts to break constraints/vsites before running the simulation

Better logging and reporting
- T* reporters
- reporter that reports the number of fragments
- reporter that reports reactions (including frame id)
- reporter that dumps T* every D/M run (including frame id)

Benchmarking polyurethane - 3 variables to create different benchmarks
- D/M frequency
- dilution
- system size
- aim for runs that finish within 5 minutes
- make pie charts for D/M, reinitialize and MD steps based on logs

Optimizations
- do not add S* forces that are never going to change (reduced memory usage when dilute systems)

Rate correctness
- limiter should not depend on the order of atoms in the topology
- reaction rates should be independent of how often the D/M algorithm is done

Behavior
- only update atom types when [rx_update] (?)

Tooling
- write itp's at the end
- checkpoints for simulations

Compatibility
- position restraints
- dihedral type 4
- pbc that are not 90 degree angle (?)
- fix pairs
- fix vsite1

Testing
- Unit test D/M algorithm
