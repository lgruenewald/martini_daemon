#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import BondReporter

sim = simulation.DaemonSimulation("system.top", "system.gro",
                                  reporters=[BondReporter])
sim.simulate()
