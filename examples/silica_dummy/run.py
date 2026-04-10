#!/usr/bin/env python3

from martini_daemon import (
    FragCountReporter,
    ReactionReporter,
    Simulation,
    ToptrajReporter,
    XTCReporter,
)

sim = Simulation(
    "system.top",
    "system.gro",
    md_steps=1000000,
    reporters=[
        ToptrajReporter(),
        XTCReporter(),
        ReactionReporter(),
        FragCountReporter(),
    ],
    dm_frequency=100,
    traj_frequency=1000,
    sim_name="out",
)
sim.get_context().minimize_energy()
sim.get_context().generate_velocities(300.0)
sim.simulate()
sim.save_geometry("out.gro")
