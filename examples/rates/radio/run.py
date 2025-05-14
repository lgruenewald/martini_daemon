#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.topstar import TopStarLogger, ReactionReporter, FragCountReporter

sim = simulation.DaemonSimulation(top_path="system.top", gro_path="system.gro", 
                                  sim_name="out",
                                  p_bar=1., friction_ps_1=2.0,
                                  reporters=[
                                      TopStarLogger(),
                                      ReactionReporter(),
                                      FragCountReporter(),
                                  ],
                                  md_steps=1000000, dm_frequency=100,
                                  xtc_frequency=10000,
                                  rate_highest_probability=0.05
                                  )
sim.simulate()
