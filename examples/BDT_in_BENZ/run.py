#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import BondReporter

sim = simulation.DaemonSimulation("system.top", "system.gro", sim_name="out", 
                                  p=None, reporters=[BondReporter],
                                  max_steps=300, steps_per_step=500)
sim.simulate()
