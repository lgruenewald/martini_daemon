#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import VMDBondReporter
from martini_daemon.reporters.topstar import TopStarLogger, ReactionReporter, FragCountReporter
from martini_daemon.reporters.checkpoint_reporter import CheckpointReporter
from martini_daemon.reporters.variables_reporter import VariablesReporter

sim = simulation.DaemonSimulation(top_path="system.top", gro_path="system.gro", 
                                  sim_name="out",
                                  p_bar=1., friction_ps_1=2.0,
                                  reporters=[
                                      VMDBondReporter(),
                                      TopStarLogger(),
                                      ReactionReporter(),
                                      FragCountReporter(),
                                      VariablesReporter(),
                                      CheckpointReporter()
                                  ],
                                  md_steps=1000000, dm_frequency=100,
                                  xtc_frequency=1000
                                  )
sim.simulate()
