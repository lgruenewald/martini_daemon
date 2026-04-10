#!/usr/bin/env python3

import openmm as mm
from martini_daemon import (
    Simulation,
    ToptrajReporter,
    VariablesReporter,
    XTCReporter,
    ReactionReporter,
    FragCountReporter,
)

sim = Simulation(
    "system.top",
    "system.gro",
    reporters=[
        ToptrajReporter(),
        VariablesReporter(),
        XTCReporter(),
        ReactionReporter(),
        FragCountReporter(),
    ],
    md_steps=10000000,
    dm_frequency=250,
    traj_frequency=500,
    platform="CUDA",
    integrator=mm.LangevinIntegrator(
        0.01 * mm.unit.picosecond, 1.0 / mm.unit.picosecond, 298.0 * mm.unit.kelvin
    ),
)
sim.context.minimize_energy()
sim.context.generate_velocities(298.0)
sim.simulate()
