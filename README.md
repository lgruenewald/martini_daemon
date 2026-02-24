# Martini Daemon

Martini Daemon is a tool facilitating template based chemical reactions in MD simulations with the [Martini force field](https://cgmartini.nl/) and the [OpenMM](https://openmm.org/) MD engine.

This is achieved by combining multiple components in one repo:
- a friendly python API for running MD simulations with reactions
- a parser for Martini `.top` files
- a thin wrapper on top of OpenMM's API facilitating bond addition and removal
- a graph matching system to find reactants
- a detection/modification algorithm to execute reaction templates

### Installation

Pre-requisite: install the right version of OpenMM with support for your GPU and install `git-lfs`.

1. Clone the repository
2. (optional) switch to the desired branch, tag or commit
3. (optional) activate the desired conda, mamba, uv or venv environment
4. Install with `pip install .` in the root directory of this repo,
   where `pyproject.toml` is located.

## Docs

Docs can be built using sphinx. Run `make html` in the `docs/` folder in the
repo. The generated docs can then be found under `docs/build/html`.

If any questions remain, feel free to open an Issue, so that we can help and also extend the documentation where necessary.

## Examples

Example systems built with martini_daemon can be found here.
The README in every example directory should provide further
information about each system.

## Tests

- `single_frame` - single point energy and force calculation tests, that verify
  that Martini is implemented correctly by comparing it to GROMACS energies.
- `graph` - tests for the graph matching algorithm
- `detection` - tests for the detection algorithm
- `modification` - tests for the modification algorithm (WIP)

# License

Martini Daemon is licensed under the Apache 2.0 license.
See LICENSE.txt for details.
