#!/usr/bin/env python

from martini_daemon.simulation import DaemonSimulation
from martini_daemon.reporters.type_reporter import TypeReporter

sim = DaemonSimulation("system.top", "system.gro", max_steps=2000,
                       steps_per_step=1000, reporters=[TypeReporter])
sim.simulate()
