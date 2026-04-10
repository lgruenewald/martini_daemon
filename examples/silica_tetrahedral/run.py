#!/usr/bin/env python3

from martini_daemon import (
    Simulation,
    ToptrajReporter,
    VariablesReporter,
    XTCReporter,
    ReactionReporter,
    FragCountReporter,
)

sim = Simulation(
    top_path="system.top",
    geom_path="system.gro",
    sim_name="out",
    reporters=[
        ToptrajReporter(),
        VariablesReporter(),
        XTCReporter(),
        ReactionReporter(),
        FragCountReporter(),
    ],
    md_steps=10000000,
    dm_frequency=250,
    traj_frequency=5000,
)
sim.context.generate_velocities(300)
sim.simulate()
