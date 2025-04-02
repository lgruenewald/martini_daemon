#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import VMDBondReporter
from martini_daemon.reporters.variables_reporter import VariablesReporter

sim = simulation.DaemonSimulation("system.top", "system.gro",
                                  reporters=[VMDBondReporter, VariablesReporter],
                                  md_steps=2000, dm_frequency=100,
                                  xtc_frequency=1000,
                                  platform="CUDA")
sim.simulate()
