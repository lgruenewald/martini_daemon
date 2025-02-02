#!/usr/bin/env python

from martini_daemon.simulation import DaemonSimulation
from martini_daemon.reporters.bond_reporter import BondReporter

sim = DaemonSimulation("chaining.top", "chaining.gro", max_steps=400,
                       steps_per_step=1000, reporters=[BondReporter])
sim.simulate()
