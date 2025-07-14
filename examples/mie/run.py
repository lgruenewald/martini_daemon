#!/usr/bin/env python3

from martini_daemon.simulation import Simulation
from martini_daemon.forces.mie import MiePotential
from martini_daemon.reporters.variables_reporter import VariablesReporter
import numpy as np
import openmm as mm

params = [(12, 6), (9, 7), (12, 4), (8, 4), (9, 6), (5, 4), (6, 3), (4, 2)]
for m, n in params:
    sim = Simulation(
        "system.top", "system.gro", md_steps=1000000, dm_frequency=0,
        traj_frequency=5000, sim_name=f"shift{m}-{n}",
        reporters=[
            VariablesReporter()
        ],
        nonbonded_force=MiePotential(n, m),
        integrator=mm.VerletIntegrator(0.02),
        coupling=[
            mm.AndersenThermostat(300, 1.0),
            mm.MonteCarloBarostat(1, 300)
        ]
    )
    sim.generate_velocities(300)
    sim.minimize_energy()
    sim.simulate()
