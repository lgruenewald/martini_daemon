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

from martini_daemon.top_parser import DaemonTopFile

etol = 1e-5  # energy relative tolerance
ftol = 1e-5  # force relative tolerance
rtol = 2e-3  # distance tolerance


def test_constraints(top, gro):
    # applies constraints and vsites and checks for position change
    platform = mm.Platform.getPlatformByName("Reference")
    system, top = DaemonTopFile(top)
    gro = mmapp.GromacsGroFile(gro)
    oldpos = gro.getPositions(True)
    system.build_context(mm.VerletIntegrator(2 * femtosecond),
                         gro.getPeriodicBoxVectors(),
                         platform=platform)
    system.set_positions(gro.getPositions(True))
    system.apply_constraints()
    newpos = system.get_state().getPositions(asNumpy=True)
    delta = np.linalg.norm(newpos - oldpos, axis=1)
    if np.any(delta > rtol):
        largest_index = np.argmax(delta)
        return False, "Constraint/VSite position deviation, largest at "\
                      f"{largest_index}: {delta[largest_index]} nm.\n"
    return True, None


def test_daemon(top, gro):
    platform = mm.Platform.getPlatformByName("Reference")
    system, top = DaemonTopFile(top)
    gro = mmapp.GromacsGroFile(gro)
    system.build_context(mm.VerletIntegrator(2 * femtosecond),
                         gro.getPeriodicBoxVectors(),
                         platform=platform)
    system.set_positions(gro.getPositions(True))
    system.apply_constraints()
    state = system.get_state()
    energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    forces = state.getForces(asNumpy=True).\
        value_in_unit(kilojoule / nanometer / mole).flatten()
    for vsite in system.vsites:
        forces[vsite * 3] = 0.
        forces[vsite * 3 + 1] = 0.
        forces[vsite * 3 + 2] = 0.
    return energy, forces


def test_martini_openmm(top, gro):
    platform = mm.Platform.getPlatformByName("Reference")
    gro = mmapp.GromacsGroFile(gro)
    box_vectors = gro.getPeriodicBoxVectors()
    top = martini.MartiniTopFile(top, periodicBoxVectors=box_vectors)
    system = top.create_system()
    integrator = mm.LangevinIntegrator(
        300, 1.0, 2 * femtosecond
    )
    sim = mmapp.Simulation(top.topology, system, integrator, platform)
    sim.context.setPositions(gro.getPositions())
    sim.context.applyConstraints(tol=1e-10)
    state = sim.context.getState(getEnergy=True, getForces=True)
    energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    forces = state.getForces(asNumpy=True).\
        value_in_unit(kilojoule / nanometer / mole).flatten()
    for vsite in top._all_vsites:
        forces[vsite * 3] = 0.
        forces[vsite * 3 + 1] = 0.
        forces[vsite * 3 + 2] = 0.
    return energy, forces


argv = sys.argv
argc = len(argv)

tests = os.listdir(os.curdir)

if argc > 1:
    tests = argv[1:]

passed = 0
failed = 0

for x in tests:
    if os.path.isdir(x):
        # Run GMX
        os.chdir(x)
        os.system("../gmxrun.sh")
        errors = []
        notes = []

        constraint_passed, msg = test_constraints("system.top", "system.gro")
        if not constraint_passed:
            errors.append(msg)
        # Run Daemon+openmm
        daemon_energy, daemon_forces = test_daemon("system.top", "system.gro")
        # run martini_openmm
        # Compare forces and energies
        with open("energy.xvg") as f:
            lines = [line for line in f]
            gmx_energy = float(lines[-1].split()[-1])

        diff = math.fabs(gmx_energy / daemon_energy - 1)

        if diff > etol:
            errors.append(
                f"Gromacs energy: {gmx_energy:.10e}\n"
                f"Daemon energy: {daemon_energy:.10e}\n"
                f"Relative difference {diff:.3e} above tolerance {etol:.2e}\n"
            )

        with open("forces.xvg") as f:
            lines = [line for line in f]
            gmx_force_line = lines[-1].split()
            gmx_forces = np.array([float(x) for x in gmx_force_line][1:])

        force_diff = np.fabs(gmx_forces - daemon_forces) / \
            (np.fabs(daemon_forces) + ftol)
        if not np.allclose(gmx_forces, daemon_forces, ftol, 0.):
            i_max = np.argmax(force_diff)
            max = force_diff[i_max]
            daemon_force = daemon_forces[i_max]
            gmx_force = gmx_forces[i_max]
            abs_diff = np.fabs(daemon_force - gmx_force)
            atom_index = i_max // 3
            atom_dim = i_max % 3
            errors.append(
                f"Large force deviation for particle {atom_index} "
                f"dimension {atom_dim}\n"
                f"absolute diff: {abs_diff:.3e}    relative diff: {max:.3e}\n"
                f"Daemon force: {daemon_force:.10e}\n"
                f"Gromacs force: {gmx_force:.10e}\n"
            )

        try:
            mo_energy, mo_forces = test_martini_openmm("system.top", "system.gro")
            diff2 = math.fabs(mo_energy / daemon_energy - 1)
            if diff2 > 1e-7:
                notes.append(
                    "Warning: martini_openmm energy different\n"
                    f"Gromacs energy: {gmx_energy:.10e}\n"
                    f"Daemon energy: {daemon_energy:.10e}\n"
                    f"Martini openmm energy: {mo_energy:.10e}.\n"
                )
            elif len(errors) > 0:
                notes.append("Martini_openmm energy matches daemon.\n")

            force_diff2 = np.fabs(mo_forces - daemon_forces) / \
                (np.fabs(daemon_forces) + 1e-7)
            if np.any(force_diff2 > 1e-7):
                i_max = np.argmax(force_diff2)
                max = force_diff2[i_max]
                daemon_force = daemon_forces[i_max]
                mo_force = mo_forces[i_max]
                abs_diff = np.fabs(daemon_force - mo_force)
                atom_index = i_max // 3
                atom_dim = i_max % 3
                notes.append(
                    f"Large force deviation for particle {atom_index} "
                    f"dimension {atom_dim}\n"
                    f"absolute diff: {abs_diff:.3e}    relative diff: {max:.3e}\n"
                    f"Daemon force: {daemon_force:.10e}\n"
                    f"Martini Openmm force: {mo_force:.10e}\n"
                )
            elif len(errors) > 0:
                notes.append("Martini openmm forces match with daemon.\n")

        except:
            notes.append("Martini openmm failed to run.\n")

        if len(errors) > 0 or len(notes) > 0:
            sys.stderr.write(f"== Test {x} ==\n")
            for error in errors:
                sys.stderr.write("\x1b[1;31m[ERROR] ")
                sys.stderr.write(error)
            for note in notes:
                sys.stderr.write("\x1b[1;33m[NOTE] ")
                sys.stderr.write(note)
            sys.stderr.write("\x1b[0m\n")

        if len(errors) > 0:
            failed += 1
        else:
            passed += 1

        # cleanup
#        os.remove("energy.xvg")
#        os.remove("forces.xvg")
        os.chdir("..")

print(f"Passed: {passed} Failed: {failed}")
