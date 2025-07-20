#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import BondReporter
from martini_daemon.reporters.topstar import ReactionReporter

sim = simulation.Simulation(
    top_path="system.top", gro_path="system.gro",
    sim_name="out",
    reporters=[
        BondReporter(),
        ReactionReporter(),
    ],
    md_steps=100000000, dm_frequency=100,
    xtc_frequency=5000,
)
sim.minimize_energy()
sim.generate_velocities(300)
sim.simulate()
