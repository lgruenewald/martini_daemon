#!/usr/bin/env python3

from martini_daemon import old_simulation
from martini_daemon.old_reporters.bond_reporter import BondReporter
from martini_daemon.old_reporters.topstar import ReactionReporter
from martini_daemon.old_components.reaction_sensitive_integrator import ReactionSensitiveLangevinIntegrator

sim = simulation.Simulation(
    top_path="system.top", gro_path="system.gro",
    sim_name="out",
    reporters=[
        BondReporter(),
        ReactionReporter(molid=True),
    ],
    md_steps=100000000, dm_frequency=100,
    xtc_frequency=5000,
    integrator=ReactionSensitiveLangevinIntegrator(0.02, 298, 1., 4, 100)
)
sim.minimize_energy()
sim.generate_velocities(300)
sim.simulate()
