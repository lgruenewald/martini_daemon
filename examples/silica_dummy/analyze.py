#!/usr/bin/env python3

from martini_daemon.old_reporters.bond_reporter import read_bonds
from martini_daemon.old_helpers.vmd import generate_vmd_readable_bonds

bonds = read_bonds("out.bonds")
generate_vmd_readable_bonds(bonds, "whole.xtc", "whole.z")
