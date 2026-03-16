#!/usr/bin/env python3

from martini_daemon.old_reporters.bond_reporter import read_bonds
from martini_daemon.old_helpers.vmd import generate_vmd_readable_bonds

n_frames, n_atoms, bonds = read_bonds("out.bonds")
print("n_frames", n_frames)
print("n_atoms", n_atoms)
generate_vmd_readable_bonds(bonds, "out.xtc", "out.z")
