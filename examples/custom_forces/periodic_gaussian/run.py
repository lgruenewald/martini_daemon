#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.variables_reporter import VariablesReporter

sim =simulation.Simulation(
    "system.top", "system.gro",
    10000000, 0, 5000, "out",
    reporters=[
        VariablesReporter()
    ]
)
sim.simulate()
