#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import VMDBondReporter
from martini_daemon.reporters.topstar import ReactionReporter

sim = simulation.DaemonSimulation(top_path="system.top", gro_path="system.gro",
                                  sim_name="out",
                                  reporters=[
                                      VMDBondReporter(),
                                      ReactionReporter(),
                                  ],
                                  md_steps=100000000, dm_frequency=100,
                                  xtc_frequency=5000,
                                  platform="CUDA"
                                  )
sim.simulate()
