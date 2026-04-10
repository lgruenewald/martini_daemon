#!/usr/bin/env python3

import openmm as mm
from openmm.unit import kelvin, picosecond  # ty: ignore[unresolved-import]

from martini_daemon import (
    FragCountReporter,
    ReactionReporter,
    Simulation,
    ToptrajReporter,
    VariablesReporter,
    XTCReporter,
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
        0.01 * picosecond,
        1.0 / picosecond,
        298.0 * kelvin,  # ty: ignore[unsupported-operator]
    ),
)
sim.get_context().minimize_energy()
sim.get_context().generate_velocities(298.0)
sim.simulate()
