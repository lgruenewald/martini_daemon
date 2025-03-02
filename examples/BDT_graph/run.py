#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import BondReporter
from martini_daemon.reporters.topstar import TopStarLogger, ReactionReporter, FragCountReporter

sim = simulation.DaemonSimulation("system.top", "system.gro", sim_name="out",
                                  p=None,
                                  reporters=[BondReporter, TopStarLogger,
                                             ReactionReporter,
                                             FragCountReporter],
                                  max_steps=1000, steps_per_step=100)
sim.simulate()
