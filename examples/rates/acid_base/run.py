#!/usr/bin/env python3

from martini_daemon import old_simulation
from martini_daemon.old_reporters.topstar import ReactionReporter, FragCountReporter

sim = simulation.Simulation(
    top_path="system.top", gro_path="system.gro",
    sim_name="out",
    reporters=[
        ReactionReporter(),
        FragCountReporter(),
    ],
    md_steps=1000000, dm_frequency=100,
    xtc_frequency=10000,
    rate_highest_probability=0.005
)
sim.simulate()
