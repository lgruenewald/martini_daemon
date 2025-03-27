#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import BondReporter
from martini_daemon.reporters.all_bond_reporter import AllBondReporter
from martini_daemon.reporters.topstar import ReactionReporter

sim = simulation.DaemonSimulation("system.top", "system.gro", sim_name="out",
                                  reporters=[BondReporter, ReactionReporter, AllBondReporter],
                                  max_steps=1000000, steps_per_step=100,
                                  friction=2.0, T_type="langevin",
                                  xtc_every=50, platform="CUDA")
sim.simulate()
