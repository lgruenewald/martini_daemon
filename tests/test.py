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

etol = 1e-4  # energy relative tolerance
ftol = 1e-4  # force relative tolerance
rtol = 2e-3  # distance tolerance


def test_constraints(top, gro):
    # applies constraints and vsites and checks for position change
    system, top = DaemonTopFile(top)
    gro = mmapp.GromacsGroFile(gro)
    oldpos = gro.getPositions(True)
    system.build_context(mm.VerletIntegrator(20 * femtosecond),
                         gro.getPeriodicBoxVectors())
    system.set_positions(gro.getPositions(True))
    system.apply_constraints()
    newpos = system.get_state().getPositions(asNumpy=True)
    delta = np.linalg.norm(newpos - oldpos, axis=1)
    if np.any(delta > rtol):
        largest_index = np.argmax(delta)
        print("Constraint/VSite position deviation, largest at "
              f"{largest_index}: {delta[largest_index]} nm.")
        return False
    return True


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


argv = sys.argv
argc = len(argv)

tests = os.listdir(os.curdir)

if argc > 1:
    tests = argv[1:]

for x in tests:
    if os.path.isdir(x):
        # Run GMX
        print(f"Test {x}:")
        os.chdir(x)
        os.system("./gmxrun.sh")
        constraint_passed = test_constraints("system.top", "system.gro")
        # Run Daemon+openmm
        daemon_energy = test_daemon("system.top", "system.gro")
        # run martini_openmm
        # Compare forces and energies
        with open("energy.xvg") as f:
            lines = [line for line in f]
            gmx_energy = float(lines[-1].split()[-1])

        diff = math.fabs(gmx_energy / daemon_energy - 1)

        gromacs_passed = False
        if diff > etol:
            print("Gromacs energy:", gmx_energy)
            print("Daemon energy:", daemon_energy)
        else:
            gromacs_passed = True

        try:
            martini_openmm_energy = test_martini_openmm("system.top", "system.gro")
            diff2 = math.fabs(martini_openmm_energy / daemon_energy - 1)
            if diff2 > etol:
                print("Warning: martini_openmm energy different")
                if gromacs_passed:
                    print("Gromacs energy:", gmx_energy)
                    print(f"Daemon energy: {daemon_energy}")
                print(f"Martini openmm energy: {martini_openmm_energy}")
        except:
            print("Note - martini openmm comparison failed")

        if not constraint_passed or not gromacs_passed:
            print("FAIL")
        else:
            print("PASS")

        # cleanup
        os.remove("energy.xvg")
        os.remove("forces.xvg")
        os.chdir("..")
    else:
        print(f"Skipping {x}, not a directory.")
