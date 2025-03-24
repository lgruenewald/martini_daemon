#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import BondReporter
from martini_daemon.reporters.topstar import ReactionReporter

sim = simulation.DaemonSimulation("system.top", "system.gro", sim_name="out",
                                  reporters=[BondReporter, ReactionReporter],
                                  max_steps=40000, steps_per_step=250,
                                  friction=1.0, T_type="langevin",
                                  xtc_every=50)
sim.simulate()
