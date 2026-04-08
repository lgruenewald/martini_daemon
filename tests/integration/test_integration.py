from martini_daemon import Simulation, XTCReporter, ToptrajReporter, LocalMinimizer, LocalGradientDescent, VariablesReporter,  ReactionReporter, FragCountReporter, ReactionEnergyReporter
import os
import glob
import pytest

@pytest.fixture
def rootdir(request):
    return os.path.dirname(request.path)

tests = [
    "silica_dummy"
]

@pytest.mark.parametrize("x", tests)
def test_integration(x, rootdir):
    os.chdir(rootdir)
    os.chdir(x)

    sim = Simulation(
        "system.top", "system.gro",
        2500,
        [
            XTCReporter(),
            ToptrajReporter(),
            LocalMinimizer(
                minimizer=LocalGradientDescent(
                    initial_step_size_nm=0.1,
                    etol=0.0001,
                    smoothing_factor=0.1
                ),
                minimization_steps=500,
                r_movable=1.,
                whole_molecule=True,
                # technically there aren't any constraints in some systems, such as silica dummy
                harmonic_constraints=True,
                report_every=10
            ),
            VariablesReporter(),
            ReactionReporter(),
            FragCountReporter(),
            ReactionEnergyReporter(write_coords=True, ext=".xyz"),
        ],
        250, 1000,
        sim_name="out"
    )
    sim.context.minimize_energy()
    sim.context.generate_velocities(298)
    sim.simulate()

    for filename in glob.glob("./out*"):
        os.remove(filename)
    for filename in glob.glob("./#*"):
        os.remove(filename)