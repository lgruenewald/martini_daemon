#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import BondReporter
from martini_daemon.reporters.atom_reporter import AtomReporter
from martini_daemon.reporters.variables_reporter import VariablesReporter
from martini_daemon.reporters.topstar import FragCountReporter, ReactionReporter
from martini_daemon.reporters.checkpoint_reporter import CheckpointReporter

sim = simulation.DaemonSimulation(
    "system.top", "system.gro",
    reporters=[
        BondReporter(),
        AtomReporter(),
        VariablesReporter(),
        ReactionReporter(),
        FragCountReporter(),
        CheckpointReporter(100000)
    ],
    md_steps=1000000, dm_frequency=100,
    xtc_frequency=500,
    platform="CUDA"
)
sim.simulate()
