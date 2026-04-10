#!/usr/bin/env python3

from martini_daemon import Simulation, ReactionReporter, FragCountReporter, XTCReporter, VariablesReporter

sim = Simulation(
    top_path="system.top", geom_path="system.gro",
    sim_name="out",
    reporters=[
        ReactionReporter(),
        FragCountReporter(),
        XTCReporter(),
        VariablesReporter(),
    ],
    md_steps=1000000, dm_frequency=100,
    traj_frequency=10000,
)
sim.simulate()
