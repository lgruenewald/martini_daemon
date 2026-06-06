"""Run short simulations to explore the Martini Daemon components in conjunction and see if any exceptions arise."""

import glob
import os

import pytest

from martini_daemon import (
    CheckpointReporter,
    FragCountReporter,
    GlobalMinimizer,
    LocalGradientDescent,
    LocalMinimizer,
    ReactionEnergyReporter,
    ReactionReporter,
    Reporter,
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

variation = [
    "global", "default", "local"
]


@pytest.mark.parametrize("x", tests)
@pytest.mark.parametrize("variant", variation)
def test_integration(x: str, variant: str, rootdir: str) -> None:
    """Run a short simulation with many reporters."""
    os.chdir(rootdir)
    os.chdir(x)

    reporters: list[Reporter] = [
        TrajectoryReporter(),
        TopTrajReporter(),
        VariablesReporter(),
        ReactionReporter(),
        FragCountReporter(),
        ReactionEnergyReporter(write_coords=True, ext=".xyz"),
        CheckpointReporter(interval=2000, n_checkpoints=5),
    ]

    match variant:
        case "local":
            reporters.append(
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
                )
            )

        case "global":
            reporters.append(
                GlobalMinimizer(
                    minimization_steps=50,
                    restraint_force_global=1000,
                    restraint_force_local=50,
                    r_movable=0.0,
                    whole_molecule=False,
                )
            )

        case "_":
            pass

    with Simulation(
        "system.top",
        "system.gro",
        5000,
        reporters,
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
