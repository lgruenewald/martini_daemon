#!/usr/bin/env python3

from martini_daemon import (
    FragCountReporter,
    ReactionReporter,
    Simulation,
    TopTrajReporter,
    TrajectoryReporter,
    FragmentsDump
)

sim = Simulation(
    "system.top",
    "system.gro",
    md_steps=1000000,
    reporters=[
        TopTrajReporter(),
        TrajectoryReporter(".xtc"),
        ReactionReporter(),
        FragCountReporter(),
        FragmentsDump()
    ],
    dm_frequency=100,
    traj_frequency=1000,
    sim_name="out",
)
sim.context.minimize_energy()
sim.context.generate_velocities(300.0)
sim.simulate()
sim.save_geometry("out.gro")
