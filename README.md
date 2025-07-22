# Martini Daemon

Readme last updated: July 21 2025

### Installation

Pre-requisite: install the right version of OpenMM with support for your GPU

1. Clone the repository
2. (optional) switch to the desired branch, tag or commit
3. (optional) activate the desired conda, mamba or venv environment
4. Install with `pip install .` in the root directory of this repo,
   where `pyproject.toml` is located.

## Docs

You can find the User guide (`user_guide.md`) and Reference
(`reference.md`) here.

## Examples

Example systems built with martini_daemon can be found here. The README in
every example directory should provide further information about each
system.

## Tests

- `single_frame` - single point energy and force calculation tests, that verify
  that Martini is implemented correctly by comparing it to GROMACS energies.
- `graph` - tests for the graph matching algorithm
- `detection` - tests for the detection algorithm
