#!/usr/bin/env python3

from martini_daemon import old_simulation
from martini_daemon.old_reporters.atom_reporter import AtomReporter
from martini_daemon.old_reporters.topstar import ReactionReporter, FragCountReporter

for rate, highest_prob in [("0", 0.1), ("1", 0.1), ("1", 0.5), ("1", 1.0)]:
    for rep in range(3):
        sim = simulation.Simulation(
            top_path="system.top", gro_path="system.gro",
            sim_name=f"out_r{rate}_p{highest_prob}_rep{rep}",
            p_bar=1., friction_ps_1=2.0,
            defines={"RATE": rate},
            reporters=[
                ReactionReporter(),
                FragCountReporter(),
                AtomReporter()
            ],
            md_steps=1000000, dm_frequency=100,
            xtc_frequency=5000,
            rate_highest_probability=highest_prob
        )
        sim.simulate()
