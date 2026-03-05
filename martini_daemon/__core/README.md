Core
====

everything here is WIP, this code is still not called

- Things in Core should only depend on library dependencies, not on other things outside of core.
- Things in Core should have only a minimal set of defaults. Defaults for most Martini simulations should be elsewhere.
- Implemented here: OpenMM abstraction layer/wrapper, System* related bookkeeping, file formats and I/O. 
- See `__init__.py` on what gets exposed.