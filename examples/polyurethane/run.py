#!/usr/bin/env python3
from martini_daemon import (
    FragCountReporter,
    ReactionReporter,
    Simulation,
    ToptrajReporter,
    TrajectoryReporter,
    VariablesReporter,
)

sim = Simulation(
    "system.top",
    "system.gro",
    reporters=[
        ReactionReporter(),
        FragCountReporter(),
        TrajectoryReporter(),
        VariablesReporter(),
        ToptrajReporter(),
    ],
    md_steps=100000000,
    dm_frequency=250,
    traj_frequency=50000,
    platform="CUDA",
    context_parameters={"DeviceIndex": "0"},
)
sim.context.minimize_energy()
sim.context.generate_velocities(300)
sim.simulate()
