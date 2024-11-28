#!/usr/bin/env python3

from daemon_top_parser import DaemonTopFile
import sys
import openmm as mm
import openmm.app as mmapp
from openmm.unit import kelvin, picosecond, femtosecond, kilojoule_per_mole, \
                        kilojoule, mole, nanometer
import numpy as np


argv = sys.argv
argc = len(sys.argv)

if argc != 3:
    print("Usage: ./main.py <top file> <gro file>")
    quit(1)

top_path = argv[1]
gro_path = argv[2]

system, top = DaemonTopFile(top_path)


gro = mmapp.GromacsGroFile(gro_path)

system.add_force(mm.MonteCarloBarostat(1.0, 300 * kelvin))

system.build_context(mm.LangevinIntegrator(300 * kelvin, 10.0 / picosecond,
                                           20 * femtosecond),
                     gro.getPeriodicBoxVectors())

pos = system.set_positions(gro.getPositions(True))

# dump forces and energies

state = system.get_state()
energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
forces = np.array(state.getForces()
                  .value_in_unit(kilojoule / nanometer / mole))

print("Forces: ", forces)
print("Energy: ", energy)
