#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import BondReporter
from martini_daemon.reporters.all_bond_reporter import AllBondReporter

sim = simulation.DaemonSimulation("system.top", "system.gro",
                                  reporters=[BondReporter, AllBondReporter],
                                  max_steps=1000, platform="CUDA")
sim.simulate()
