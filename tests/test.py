#!/usr/bin/env python3

import os
import sys
import openmm as mm
import openmm.app as mmapp
from openmm.unit import kelvin, picosecond, femtosecond,kilojoule_per_mole, \
                        kilojoule, mole, nanometer
import numpy as np
import math

sys.path.append("../src")
from daemon_top_parser import DaemonTopFile

atol = 1e-3
rtol = 1e-4

def test_openmm(top, gro):
    system, top = DaemonTopFile(top)
    gro = mmapp.GromacsGroFile(gro)
    system.build_context(mm.VerletIntegrator(20 * femtosecond),
                         gro.getPeriodicBoxVectors())
    system.set_positions(gro.getPositions(True))
    state = system.get_state()
    energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    return energy
    

for x in os.listdir(os.curdir):
    if os.path.isdir(x):
        # Run GMX
        os.chdir(x)
        os.system("./gmxrun.sh")
        # Run Daemon+openmm
        openmm_energy = test_openmm("system.top", "system.gro")

        # Compare forces and energies
        with open("energy.xvg") as f:
            lines = [line for line in f]
            gmx_energy = float(lines[-1].split()[-1])

        if math.fabs(gmx_energy - openmm_energy) > atol or \
                math.fabs(openmm_energy / gmx_energy - 1) > rtol:
            print(f"Failed test {x}")
            print("Gromacs energy:", gmx_energy)
            print("Openmm energy:", openmm_energy)
        else:
            print(f"Passed test {x}")

        # cleanup
        os.remove("energy.xvg")
        os.remove("forces.xvg")
        os.chdir("..")
