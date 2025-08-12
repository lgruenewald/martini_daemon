# Martini Daemon

Martini Daemon is a tool facilitating template based chemical reactions in MD simulations with the [Martini force field](https://cgmartini.nl/) and the [OpenMM](https://openmm.org/) MD engine.

### Installation

Pre-requisite: install the right version of OpenMM with support for your GPU and install `git-lfs`.

1. Clone the repository
2. (optional) switch to the desired branch, tag or commit
3. (optional) activate the desired conda, mamba or venv environment
4. Install with `pip install .` in the root directory of this repo,
   where `pyproject.toml` is located.

## Docs

Currently, there is a User Guide and a Reference Guide in the wiki section of this Github repository.

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
