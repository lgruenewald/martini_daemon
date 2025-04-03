#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import VMDBondReporter
from martini_daemon.reporters.topstar import TopStarLogger, ReactionReporter, FragCountReporter

sim = simulation.DaemonSimulation(top_path="system.top", gro_path="system.gro",
                                  sim_name="out",
                                  p_bar=1.0, friction_ps_1=2.0,
                                  reporters=[
                                      VMDBondReporter(),
                                      TopStarLogger(),
                                      ReactionReporter(),
                                      FragCountReporter()
                                  ],
                                  md_steps=1000000, dm_freq=100,
                                  xtc_every=1000,
                                  minimize_energy=False)
sim.simulate()
