#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import BondReporter
from martini_daemon.reporters.all_bond_reporter import AllBondReporter
from martini_daemon.reporters.topstar import ReactionReporter

sim = simulation.DaemonSimulation(top_path="system.top", gro_path="system.gro",
                                  sim_name="out",
                                  reporters=[
                                      BondReporter(),
                                      ReactionReporter(),
                                      AllBondReporter()
                                  ],
                                  md_steps=100000000, dm_freq=100,
                                  xtc_freq=5000,
                                  friction=2.0,
                                  platform="CUDA"
                                  )
sim.simulate()
