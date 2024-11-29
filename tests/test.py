#!/usr/bin/env python3

import os
import sys
import openmm as mm
import openmm.app as mmapp
from openmm.unit import kelvin, picosecond, femtosecond,kilojoule_per_mole, \
                        kilojoule, mole, nanometer
import numpy as np
import math
import martini_openmm as martini

sys.path.append("../src")
from daemon_top_parser import DaemonTopFile

tol = 1e-4


def test_daemon(top, gro):
    system, top = DaemonTopFile(top)
    gro = mmapp.GromacsGroFile(gro)
    system.build_context(mm.VerletIntegrator(20 * femtosecond),
                         gro.getPeriodicBoxVectors())
    system.set_positions(gro.getPositions(True))
    state = system.get_state()
    energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    return energy


def test_martini_openmm(top, gro):
    gro = mmapp.GromacsGroFile(gro)
    box_vectors = gro.getPeriodicBoxVectors()
    top = martini.MartiniTopFile(top, periodicBoxVectors=box_vectors)
    system = top.create_system()
    integrator = mm.LangevinIntegrator(
        300, 1.0, 2
    )
    sim = mmapp.Simulation(top.topology, system, integrator)
    sim.context.setPositions(gro.getPositions())
    state = sim.context.getState(getEnergy=True, getForces=True)
    return state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)


for x in os.listdir(os.curdir):
    if os.path.isdir(x):
        # Run GMX
        os.chdir(x)
        os.system("./gmxrun.sh")
        # Run Daemon+openmm
        daemon_energy = test_daemon("system.top", "system.gro")
        # run martini_openmm
        martini_openmm_energy = test_martini_openmm("system.top", "system.gro")

        # Compare forces and energies
        with open("energy.xvg") as f:
            lines = [line for line in f]
            gmx_energy = float(lines[-1].split()[-1])

        diff = math.fabs(gmx_energy / daemon_energy - 1)
        if diff > tol:
            print(f"Failed test {x}")
            print("Gromacs energy:", gmx_energy)
            print("Daemon energy:", daemon_energy)
            print("(for context)martini_openmm energy:", martini_openmm_energy)
        else:
            print(f"Passed test {x}")

        # cleanup
        os.remove("energy.xvg")
        os.remove("forces.xvg")
        os.chdir("..")
