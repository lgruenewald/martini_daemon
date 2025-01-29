#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import BondReporter

sim = simulation.DaemonSimulation("system.top", "system.gro", p=None, reporters=[BondReporter], max_steps=1000, top_logpath="top.log")
sim.simulate()
