#!/usr/bin/env python3

from martini_daemon.reporters.bond_reporter import read_bonds
from martini_daemon.helpers.vmd import generate_vmd_readable_bonds

bonds = read_bonds("out.bonds")
generate_vmd_readable_bonds(bonds, "out.xtc", "out.z")
