#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.type_reporter import TypeReporter

sim = simulation.DaemonSimulation("system.top", "system.gro", sim_name="out", 
                                  reporters=[TypeReporter],
                                  max_steps=20, steps_per_step=0)
sim.simulate()
