#!/usr/bin/env python3

from martini_daemon import old_simulation
from martini_daemon.old_reporters.bond_reporter import BondReporter
from martini_daemon.old_reporters.atom_reporter import AtomReporter
from martini_daemon.old_reporters.variables_reporter import VariablesReporter
from martini_daemon.old_reporters.topstar import FragCountReporter, ReactionReporter
from martini_daemon.old_reporters.checkpoint_reporter import CheckpointReporter

sim = simulation.Simulation(
    "system.top", "system.gro",
    reporters=[
        BondReporter(),
        AtomReporter(),
        VariablesReporter(),
        ReactionReporter(),
        FragCountReporter(),
        CheckpointReporter(100000)
    ],
    md_steps=10000000, dm_frequency=250,
    xtc_frequency=500,
    platform="CUDA",
    dt_ps=0.01
)
sim.minimize_energy()
sim.generate_velocities(300)
sim.simulate()
