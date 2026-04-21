#!/usr/bin/env python3

from martini_daemon import (
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
    ],
    md_steps=100000000,
    dm_frequency=100,
    traj_frequency=5000,
)
sim.context.minimize_energy()
sim.context.generate_velocities(300)
sim.simulate()
