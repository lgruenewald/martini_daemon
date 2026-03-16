# Martini Daemon

Martini Daemon is a tool facilitating template based chemical reactions in MD simulations with the [Martini force field](https://cgmartini.nl/) and the [OpenMM](https://openmm.org/) MD engine.

This is achieved by combining multiple components in one repo / one python package:
- a friendly Python API for running MD simulations with reactions (`simulation.py`)
- a parser for GROMACS `.top` files targeting OpenMM (`__parser`, `__forces` and `__vsites`)
- a thin wrapper on top of OpenMM's API facilitating bond addition and removal (`__core`)
- a graph matching system to find reactants (`__topstar`)
- a detection/modification algorithm to execute reaction templates (`__topstar`)

### Installation

Pre-requisites:
- Make and activate a Python virtual environment or conda environment. Python 3.13 is recommended.
- It's recommended to explicitly install the right version of OpenMM with support for your GPU.
   - e.g. `pip install openmm[cuda12]`
- Install `git-lfs` (`sudo apt install git-lfs` on Ubuntu).
- Install Cargo and Rust. Installation via [rustup](https://rustup.rs/) is recommended.

Installation process:

- Download the source code. Clone the repository and switch to the desired branch, tag or commit.
- Install with `pip install .` in the root directory of this repo, where `pyproject.toml` is located.
- (optional) Install with `pip install .[all]` if you want to run tests, benchmarks or build documentation yourself.

## Docs

Docs can be built using sphinx. Run `make html` in the `docs/` folder in the
repo. The generated docs can then be found under `docs/build/html`.

If any questions remain, feel free to open an Issue, so that we can help and also extend the documentation where necessary.

## Examples

Example systems built with martini_daemon can be found here.
The README in every example directory should provide further
information about each system.

## Tests

- `test_pbc.py` - Test the PBC handling utility functions.
- `parser` - Tests the .top parser basics. 
- `single_frame` - Single point energy and force calculation tests, that verify that Martini is implemented correctly by comparing it to GROMACS energies.
- `graph` - Tests for the graph matching algorithm.
- `detection` - Tests for the detection algorithm.
- `modification` - Tests for the modification algorithm.

# License

Martini Daemon is licensed under the Apache 2.0 license.
See LICENSE.txt for details.
