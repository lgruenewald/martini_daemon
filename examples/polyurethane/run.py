#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import VMDBondReporter
from martini_daemon.reporters.all_bond_reporter import AllBondReporter

sim = simulation.DaemonSimulation("system.top", "system.gro",
                                  reporters=[VMDBondReporter,AllBondReporter],
                                  max_steps=1000, platform="CUDA",
                                  friction=10.0)
sim.simulate()
