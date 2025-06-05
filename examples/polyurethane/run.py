#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import BondReporter
from martini_daemon.reporters.variables_reporter import VariablesReporter
from martini_daemon.reporters.topstar import FragCountReporter, ReactionReporter
from martini_daemon.reporters.checkpoint_reporter import CheckpointReporter
import freud.parallel


gpu = 0
n_threads = 10

freud.parallel.set_num_threads(n_threads)

sim = simulation.Simulation(
    "system.top", "system.gro",
    reporters=[
        BondReporter(),
        VariablesReporter(),
        ReactionReporter(),
        FragCountReporter(),
        CheckpointReporter(100000)
    ],
    md_steps=100000000, dm_frequency=250,
    xtc_frequency=50000,
    platform="CUDA",
    context_parameters={"DeviceIndex": f"{gpu}"}
)
sim.simulate()
