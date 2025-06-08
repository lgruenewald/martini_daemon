#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.variables_reporter import VariablesReporter
from martini_daemon.reporters.bond_reporter import BondReporter

sim = simulation.Simulation(
    top_path="system.top", gro_path="system.gro",
    sim_name="out",
    reporters=[
        VariablesReporter()
    ],
    md_steps=1000000, dm_frequency=100,
    xtc_frequency=1000,
)
sim.generate_velocities(300)
sim.simulate()
