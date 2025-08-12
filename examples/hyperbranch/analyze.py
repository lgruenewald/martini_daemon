#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import read_bonds
from martini_daemon.helpers.vmd import generate_vmd_readable_bonds

# just to get the initial molecules from the .top file
init_mol = simulation.Simulation(
    "system.top", "system.gro",
    # it still writes logs by default
    sim_name="#tmp",
).top.initial_molecules

print("loading bonds data")
bonds_data = read_bonds("out.bonds")
generate_vmd_readable_bonds(bonds_data, "out.xtc", "out.z")
