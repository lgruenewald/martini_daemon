#!/usr/bin/env python3

from martini_daemon import (
    ReactionReporter,
    Simulation,
    ToptrajReporter,
    VariablesReporter,
    XTCReporter,
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
    ],
    md_steps=100000000,
    dm_frequency=100,
    traj_frequency=5000,
)
sim.get_context().minimize_energy()
sim.get_context().generate_velocities(300)
sim.simulate()
