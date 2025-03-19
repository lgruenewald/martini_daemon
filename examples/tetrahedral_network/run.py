#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import BondReporter
from martini_daemon.reporters.topstar import TopStarLogger, ReactionReporter, FragCountReporter

sim = simulation.DaemonSimulation("system.top", "system.gro", sim_name="out",
                                  p=1., friction=2.0,
                                  reporters=[BondReporter, TopStarLogger,
                                             ReactionReporter,
                                             FragCountReporter],
                                  max_steps=2500, steps_per_step=5000)
sim.simulate()
