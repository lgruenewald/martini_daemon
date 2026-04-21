#!/usr/bin/env python3

import openmm as mm

from martini_daemon import (
    FragCountReporter,
    ReactionReporter,
    Simulation,
    TopTrajReporter,
    TrajectoryReporter,
    VariablesReporter,
)

sim = Simulation(
    "system.top",
    "system.gro",
    reporters=[
        TopTrajReporter(),
        VariablesReporter(),
        TrajectoryReporter(),
        ReactionReporter(),
        FragCountReporter(),
    ],
    md_steps=10000000,
    dm_frequency=250,
    traj_frequency=500,
    platform="CUDA",
    integrator=mm.LangevinIntegrator(
        0.01,
        1.0,
        298.0,
    ),
)
sim.context.minimize_energy()
sim.context.generate_velocities(298.0)
sim.simulate()
