#!/usr/bin/env python3
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
        ReactionReporter(),
        FragCountReporter(),
        XTCReporter(),
        VariablesReporter(),
        ToptrajReporter(),
    ],
    md_steps=100000000,
    dm_frequency=250,
    traj_frequency=50000,
    platform="CUDA",
    context_parameters={"DeviceIndex": "0"},
)
sim.get_context().minimize_energy()
sim.get_context().generate_velocities(300)
sim.simulate()
