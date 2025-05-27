#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import BondReporter
from martini_daemon.reporters.variables_reporter import VariablesReporter
from martini_daemon.reporters.topstar import FragCountReporter, ReactionReporter
from martini_daemon.reporters.atom_reporter import AtomReporter

sim = simulation.DaemonSimulation(
    "system.top", "system.gro",
    reporters=[
        BondReporter(),
        VariablesReporter(),
        ReactionReporter(),
        FragCountReporter(),
        AtomReporter()
    ],
    md_steps=2000000, dm_frequency=100,
    xtc_frequency=5000,
    platform="CUDA"
)
sim.simulate()
