#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import BondReporter
from martini_daemon.reporters.topstar import ReactionReporter

sim = simulation.DaemonSimulation("system.top", "system.gro", sim_name="out",
                                  reporters=[BondReporter, ReactionReporter],
                                  max_steps=50000, steps_per_step=100,
                                  langevin_friction=1.0)
sim.simulate()
