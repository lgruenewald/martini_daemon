#!/usr/bin/env python3

from martini_daemon import old_simulation
from martini_daemon.old_reporters.bond_reporter2 import BondReporter
from martini_daemon.old_reporters.topstar import ReactionReporter, FragCountReporter
from martini_daemon.old_reporters.checkpoint_reporter import CheckpointReporter
from martini_daemon.old_reporters.variables_reporter import VariablesReporter

sim = simulation.Simulation(
    top_path="system.top", gro_path="system.gro", 
    sim_name="out",
    reporters=[
        BondReporter(),
        ReactionReporter(),
        FragCountReporter(),
        VariablesReporter(),
        CheckpointReporter(100000)
    ],
    md_steps=10000000, dm_frequency=250,
    xtc_frequency=5000
)
sim.generate_velocities(300)
sim.simulate()
