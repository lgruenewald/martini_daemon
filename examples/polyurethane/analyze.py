#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.helpers.monomer import read_monomer_graph, generate_mapping, loop_analysis
from martini_daemon.reporters.bond_reporter import read_bonds
from martini_daemon.helpers.vmd import generate_vmd_readable_bonds

# just to get the initial molecules from the .top file
init_mol = simulation.DaemonSimulation(
    "system.top", "system.gro",
    # it still writes logs by default
    sim_name="#tmp",
    minimize_energy=False
).top.initial_molecules

print("loading bonds data")
bonds_data = read_bonds("out.bonds")
print("generating mapping")
n_monomer, n_atom, mapping = generate_mapping(init_mol)
print("generating monomer data")
monomer_data = read_monomer_graph([bonds_data[-1]], mapping)
print("finding loops")
loops = loop_analysis(n_monomer, monomer_data, max_loop=6, progress=True)
print("Last frame bonds and loops:")
print(monomer_data[-1])
print(loops[-1])

generate_vmd_readable_bonds(bonds_data, "out.xtc", "out.z")
