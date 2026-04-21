#!/usr/bin/env python3

from martini_daemon import (
    FragCountReporter,
    ReactionReporter,
    Simulation,
    TopTrajReporter,
    TrajectoryReporter,
    VariablesReporter,
)

sim = Simulation(
    top_path="system.top",
    geometry="system.gro",
    sim_name="out",
    reporters=[
        TopTrajReporter(),
        VariablesReporter(),
        TrajectoryReporter(),
        ReactionReporter(),
        FragCountReporter(),
    ],
    md_steps=10000000,
    dm_frequency=250,
    traj_frequency=5000,
)
sim.context.generate_velocities(300)
sim.simulate()
