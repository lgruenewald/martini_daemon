#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.atom_reporter import AtomReporter
from martini_daemon.reporters.topstar import ReactionReporter, FragCountReporter

sim = simulation.DaemonSimulation(
    top_path="system.top", gro_path="system.gro",
    sim_name="out",
    p_bar=1., friction_ps_1=2.0,
    reporters=[
        ReactionReporter(),
        FragCountReporter(),
        AtomReporter()
    ],
    md_steps=100000, dm_frequency=100,
    xtc_frequency=500,
    rate_highest_probability=0.1
)
sim.simulate()
