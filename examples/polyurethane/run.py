#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import VMDBondReporter
from martini_daemon.reporters.all_bond_reporter import AllBondReporter
from martini_daemon.reporters.topstar import FragCountReporter

sim = simulation.DaemonSimulation("system.top", "system.gro",
                                  reporters=[VMDBondReporter],
                                  max_steps=2000, platform="CUDA",
                                  steps_per_step=100,
                                  friction=1.0)
sim.simulate()
