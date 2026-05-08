# Martini Daemon

Martini Daemon is a tool facilitating template based chemical reactions in MD simulations with the [Martini force field](https://cgmartini.nl/) and the [OpenMM](https://openmm.org/) MD engine.

This is achieved by combining multiple components in one repo / one python package:
- Thin wrapper on top of OpenMM's API (`__core`)
- Parser for GROMACS `.top` files (`__parser`, `__forces` and `__vsites`)
- Friendly Python API for running MD simulations with reactions (`simulation.py`)
- User experience somewhere between GROMACS and OpenMM, to create a familiar workflow for running any Martini simulation with OpenMM.
- Graph matching and detection/modification algorithms to facilitate template based reactions (`__topstar`).

## Installation (from source)

- Install `git-lfs` (`sudo apt install git-lfs` on Ubuntu).
- Install Cargo and Rust. Installation via [rustup](https://rustup.rs/) is recommended.
- Clone the repository and switch to the desired branch, tag or commit. Enter the directory.
- Make and activate a Python virtual environment. Python 3.13 is recommended.
- It's recommended to explicitly [install the right version of OpenMM](https://docs.openmm.org/latest/userguide/application/01_getting_started.html#installing-openmm) with support for your GPU.
   - e.g. `pip install openmm[cuda12]` or `pip install openmm[hip7]`
   - (optional) after installing, verify which platforms are available with `python -m openmm.testInstallation`
- pip: Install using `pip install -e .`
- uv: Install maturin using `uv tool install maturin` and then run `maturin develop -r`

# Optional dependencies

## Additional trajectory formats

By default, Martini Daemon uses [molly](https://github.com/ma3ke/molly) as its XTC reader and writer.
Optional dependencies can enable other trajectory backends.
Here is a list of optional dependency tags, based on what's currently possible:

- [trr] - Gromacs .trr files, uses [mdtraj](https://www.mdtraj.org). Full precision, stores velocities.

## Docs

Docs can be built using sphinx. Install additional dependencies with `pip install .[docs]` first.

Run `make html` in the `docs/` folder in the
repo. The generated docs can then be found under `docs/build/html`.

If any questions remain, feel free to open an Issue, so that we can help and also extend the documentation where necessary.

## Tests

To run the test suite:

```
# Gromacs in double precision is required, this should be installed first
source /usr/local/gromacs-2026.1-double/bin/GMXRC
# (in the root of the repo)
# make sure additional test dependencies are installed
pip install .[test]
# make sure everything is recompiled
maturin develop
# run the python tests
pytest .
```

The python tests can be found in the `tests` folder in the repo, containing the following types of tests:
- `parser` - Tests the .top parser basics.
- `single_frame` - Single point energy and force calculation tests,
  that verify that Martini is implemented correctly by comparing it to GROMACS energies.
  Inspired by the setup that validates [martini_openmm](https://github.com/maccallumlab/martini_openmm) against GROMACS.
- `graph` - Tests for the graph matching algorithm.
- `detection` - Tests for the detection algorithm.
- `modification` - Tests for the modification algorithm.
- `integration` - Runs a short reactive simulation with various reporters.
  Does not automatically verify output at the moment, doing that is the job of the other tests.
- `formats` - Runs tests for geometry, trajectory and topology trajectory formats.
- `replay` - Test checkpoint loading, replaying reactions and comparing subsequent energies against GROMACS.
- `test_periodic_box` - Fuzz the PeriodicBox implementation against mdtraj.
- `test_bond_graph` - Test the class `BondGraph`, which features things, such as pbc whole.

Note: Running all the tests may take ~10-15 minutes.

If a test does not pass, please open an Issue.

Some single_frame tests have looser tolerances, this is documented at the top of `tests/single_frame/test_single_frame.py`.

## Development

Optional development tools are installed using `pip install .[dev]`. Here
is an overview of which these are:

- `pylsp` - Python LSP server for autocomplete.
- `ruff` - the linter and formatter used. Run using `ruff check` and `ruff format`.
- `ty` - type checker. Run using `ty check`.

If these dependencies are installed,
linting, type checking and formatting can be done with a single command:

```
ty check && ruff check --fix && ruff format
```

To re-generate the type stubs for the rust parts, run `cargo run --bin stub_gen`.
This only needs to be done when changing the rust part of the code, as the .pyi
file (`martini_daemon/__rust/__init__.pyi`) is commited to the repository.

# License

Martini Daemon is licensed under the Apache 2.0 license.
See `LICENSE.txt` for details.
