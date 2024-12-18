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
    platform = mm.Platform.getPlatformByName("Reference")
    system, top = DaemonTopFile(top)
    gro = mmapp.GromacsGroFile(gro)
    oldpos = gro.getPositions(True)
    system.build_context(mm.VerletIntegrator(20 * femtosecond),
                         gro.getPeriodicBoxVectors(),
                         platform=platform)
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
    platform = mm.Platform.getPlatformByName("Reference")
    system, top = DaemonTopFile(top)
    gro = mmapp.GromacsGroFile(gro)
    system.build_context(mm.VerletIntegrator(20 * femtosecond),
                         gro.getPeriodicBoxVectors(),
                         platform=platform)
    system.set_positions(gro.getPositions(True))
    state = system.get_state()
    energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    forces = state.getForces(asNumpy=True).\
        value_in_unit(kilojoule / nanometer / mole).flatten()
    return energy, forces


def test_martini_openmm(top, gro):
    platform = mm.Platform.getPlatformByName("Reference")
    gro = mmapp.GromacsGroFile(gro)
    box_vectors = gro.getPeriodicBoxVectors()
    top = martini.MartiniTopFile(top, periodicBoxVectors=box_vectors)
    system = top.create_system()
    integrator = mm.LangevinIntegrator(
        300, 1.0, 2
    )
    sim = mmapp.Simulation(top.topology, system, integrator, platform)
    sim.context.setPositions(gro.getPositions())
    state = sim.context.getState(getEnergy=True, getForces=True)
    energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    forces = state.getForces(asNumpy=True).\
        value_in_unit(kilojoule / nanometer / mole).flatten()
    return energy, forces


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
        os.system("../gmxrun.sh")
        constraint_passed = test_constraints("system.top", "system.gro")
        # Run Daemon+openmm
        daemon_energy, daemon_forces = test_daemon("system.top", "system.gro")
        # run martini_openmm
        # Compare forces and energies
        with open("energy.xvg") as f:
            lines = [line for line in f]
            gmx_energy = float(lines[-1].split()[-1])

        diff = math.fabs(gmx_energy / daemon_energy - 1)

        gromacs_passed = True
        if diff > etol:
            print(f"Gromacs energy: {gmx_energy:.10e}")
            print(f"Daemon energy: {daemon_energy:.10e}")
            gromacs_passed = False

        with open("forces.xvg") as f:
            lines = [line for line in f]
            gmx_force_line = lines[-1].split()
            gmx_forces = np.array([float(x) for x in gmx_force_line][1:])

        force_passed = True
        force_diff = np.fabs(gmx_forces / daemon_forces - 1.)
        if np.any(force_diff > ftol):
            print("Unmatched forces")
            i_max = np.argmax(force_diff)
            max = force_diff[i_max]
            daemon_force = daemon_forces[i_max]
            gmx_force = gmx_forces[i_max]
            abs_diff = np.fabs(daemon_force - gmx_force)
            atom_index = i_max // 3
            atom_dim = i_max % 3
            print(f"Largest deviation for particle {atom_index} "
                  f"dimension {atom_dim}: {max}")
            print(f"Daemon force: {daemon_force:.10e}")
            print(f"Gromacs force: {gmx_force:.10e}")
            print(f"Absolute difference: {abs_diff:.3e}")
            force_passed = False

        try:
            mo_energy, mo_forces = test_martini_openmm("system.top", "system.gro")
            diff2 = math.fabs(mo_energy / daemon_energy - 1)
            if diff2 > 1e-7:
                print("Warning: martini_openmm energy different")
                if gromacs_passed:
                    print(f"Gromacs energy: {gmx_energy:.10e}")
                    print(f"Daemon energy: {daemon_energy:.10e}")
                print(f"Martini openmm energy: {mo_energy:.10e}")
            elif not gromacs_passed:
                print("Note - martini_openmm energy same as daemon.")

            force_diff2 = np.fabs(mo_forces / daemon_forces - 1.)
            if np.any(force_diff2 > 1e-7):
                print("Warning - Unmatched martini_openmm force")
                i_max = np.argmax(force_diff2)
                max = force_diff2[i_max]
                daemon_force = daemon_forces[i_max]
                mo_force = mo_forces[i_max]
                abs_diff = np.fabs(daemon_force - mo_force)
                atom_index = i_max // 3
                atom_dim = i_max % 3
                print(f"Largest deviation for particle {atom_index} "
                      f"dimension {atom_dim}: {max}")
                print(f"Daemon force: {daemon_force:.10e}")
                print(f"Martini Openmm force: {mo_force:.10e}")
                print(f"Absolute difference: {abs_diff:.3e}")
            elif not force_passed:
                print("Note - martini openmm forces match with daemon.")

        except:
            print("Note - martini openmm comparison failed")

        if not constraint_passed or not gromacs_passed or not force_passed:
            print("FAIL")
        else:
            print("PASS")

        # cleanup
        os.remove("energy.xvg")
        os.remove("forces.xvg")
        os.chdir("..")
    else:
        print(f"Skipping {x}, not a directory.")
