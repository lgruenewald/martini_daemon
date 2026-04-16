#!/usr/bin/env python3

from martini_daemon import (
    FragCountReporter,
    ReactionReporter,
    Simulation,
    ToptrajReporter,
    TrajectoryReporter,
    VariablesReporter,
)

for rep in range(3):
    sim = Simulation(
        top_path="system.top",
        geometry="system.gro",
        sim_name=f"out_rep{rep}",
        defines={"RATE": "1"},
        reporters=[
            ReactionReporter(),
            FragCountReporter(),
            TrajectoryReporter(),
            VariablesReporter(),
            ToptrajReporter(),
        ],
        md_steps=1000000,
        dm_frequency=100,
        traj_frequency=5000,
    )
    sim.simulate()
