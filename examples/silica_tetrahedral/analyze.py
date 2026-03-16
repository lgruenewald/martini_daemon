#!/usr/bin/env python3

from martini_daemon.old_reporters.bond_reporter2 import read_bonds
from martini_daemon.old_helpers.vmd import generate_vmd_readable_bonds

n_atoms, bonds = read_bonds("out.bonds")
generate_vmd_readable_bonds(bonds, "out.xtc", "out.z")
