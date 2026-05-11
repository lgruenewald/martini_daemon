"""Run short simulations to explore the Martini Daemon components in conjunction and see if any exceptions arise."""

import glob
import os

import pytest

from martini_daemon import (
    CheckpointReporter,
    FragCountReporter,
    LocalGradientDescent,
    LocalMinimizer,
    ReactionEnergyReporter,
    ReactionReporter,
    Simulation,
    TopTrajReporter,
    TrajectoryReporter,
    VariablesReporter,
)


@pytest.fixture
def rootdir(request: pytest.FixtureRequest) -> str:
    """Return the root directory for this test."""
    return os.path.dirname(request.path)


tests = ["silica_dummy"]


@pytest.mark.parametrize("x", tests)
def test_integration(x: str, rootdir: str) -> None:
    """Run a short simulation with many reporters."""
    os.chdir(rootdir)
    os.chdir(x)

    with Simulation(
        "system.top",
        "system.gro",
        5000,
        [
            TrajectoryReporter(),
            TopTrajReporter(),
            LocalMinimizer(
                minimizer=LocalGradientDescent(
                    initial_step_size_nm=0.1, etol=0.0001, smoothing_factor=0.1
                ),
                minimization_steps=500,
                r_movable=1.0,
                whole_molecule=True,
                # technically there aren't any constraints in some systems, such as silica dummy
                harmonic_constraints=True,
                report_every=10,
            ),
            VariablesReporter(),
            ReactionReporter(),
            FragCountReporter(),
            ReactionEnergyReporter(write_coords=True, ext=".xyz"),
            CheckpointReporter(interval=2000, n_checkpoints=5),
        ],
        250,
        1000,
        sim_name="out",
    ) as sim:
        sim.context.minimize_energy()
        sim.context.generate_velocities(298)
        sim.simulate()

    for filename in glob.glob("./out*"):
        os.remove(filename)
    for filename in glob.glob("./#*"):
        os.remove(filename)
