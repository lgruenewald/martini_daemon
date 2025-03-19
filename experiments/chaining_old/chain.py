# OpenMM example of a liquid of n particles, where harmonic bonds can be formed
# if two beads come within a cutoff distance

import openmm as mm
from openmm import Vec3
from openmm.unit import kelvin, picosecond, nanometer
from sys import stderr, stdout
import random
import math
import numpy as np
import numba as nb
import time

# input config, TODO unhardcode, unglobal
steps_per_frame = 1000  # do frames (= bond calculations) every n steps
mass = 10.0
charge = 0.0
sigma = 1.  # nm
epsilon = 1.  # kJ
temperature = 300 * kelvin
friction = 10 / picosecond
step_size = 0.005 * picosecond
size = 10.  # nm, period box
end_time = 5000.0 * picosecond
bond_length = 0.25  # nm
bond_force = 1000.0
bond_cutoff = 0.35  # nm, below this threshold a bond is formed
n_particles = 250

bond_angle = 2 / 3 * math.pi
angle_force = 1000.0


# Interface with openmm - moved here so bond and angle forces, integrators can
# be swapped out easier
def gen_force():
    f = mm.HarmonicBondForce()
    f.setUsesPeriodicBoundaryConditions(True)
    return f


def add_force(force, id0, id1):
    force.addBond(id0, id1, bond_length, bond_force)


def gen_angle_force():
    f = mm.HarmonicAngleForce()
    f.setUsesPeriodicBoundaryConditions(True)
    return f


def add_angle(force, id0, id1, id2):
    force.addAngle(id1, id0, id2, bond_angle, angle_force)


def get_integrator():
    return mm.LangevinMiddleIntegrator(temperature, friction, step_size)
#    return VerletIntegrator(step_size)


def make_bond(bonds, angles, active, i0, i1):
    # Changes T* to create a bond between i0 and i1, updates forces
    add_force(bonds, i0, i1)

    # Creating a new chain, disabling one of the ends so it forms long chains
    if active[i0] == -1:
        active[i0] = i1
    else:
        add_angle(angles, i0, active[i0], i1)
        active[i0] = -2
    if active[i1] == -1:
        active[i1] = i0
    else:
        add_angle(angles, i1, active[i1], i0)
        active[i1] = -2


@nb.jit
def pdist(v1, v2, size):
    # periodic boundary adjusted difference between two Vec3
    # based on ReferenceForce.cpp getDeltaRPeriodic from OpenMM

    diff = v1 - v2
    base = np.floor(diff / size + 0.5) * size
    return np.sqrt(np.dot(diff-base, diff-base))


assert (np.abs(pdist(np.array([0., 0., 1.]), np.array([0., 0., 4.]), 5.) -
               2.) < 0.00001)
assert (np.abs(pdist(np.array([0., 0., 1.]), np.array([1., 0., 0.]), 5.) -
               np.sqrt(2)) < 0.00001)


@nb.jit
def detection_algorithm(active, pos):
    # returns a list of pair of bonds to be formed
    # limitation: only one reaction per active atom per step
    n_atoms = len(active)
    # if mask is set to 1, no longer considered for reactions in this step
    mask = np.repeat(0, n_atoms)
    to_be_formed = []
    for i in range(n_atoms):
        if active[i] == -2 or mask[i] == 1:
            continue
        neighbor = active[i]
        # only up to i, to prevent double count and self reactions
        for j in range(i):
            # skip fully bonded or ones already bonding in this step
            if active[j] == -2 or mask[j] == 1:
                continue
            new_neighbor = active[j]
            # no neighbor reactions, no neighbors-neighbor
            if neighbor == j or (neighbor >= 0 and neighbor == new_neighbor):
                continue
            dist = pdist(pos[i], pos[j], size)
            if dist < bond_cutoff:
                to_be_formed.append((i, j))
                mask[i] = 1
                mask[j] = 1
                break

    return to_be_formed


def detection_modification(bonds, angles, active, pos):
    # returns whether there was anything changes
    to_be_formed = detection_algorithm(active, pos)
    for (i, j) in to_be_formed:
        make_bond(bonds, angles, active, i, j)
        return True
    return False


def simulate():
    # context pre-setup, except particles
    system = mm.System()
    nonbond = mm.NonbondedForce()
    system.addForce(nonbond)
    bonds = gen_force()
    system.addForce(bonds)
    angles = gen_angle_force()
    system.addForce(angles)
    integrator = get_integrator()

    # particles setup
    for i in range(n_particles):
        system.addParticle(mass)
        nonbond.addParticle(charge, sigma, epsilon)

    # context setup
    context = mm.Context(system, integrator)
    context.setPeriodicBoxVectors(
        Vec3(size, 0.0, 0.0), Vec3(0.0, size, 0.0), Vec3(0.0, 0.0, size)
    )
    context.setPositions(np.random.rand(n_particles, 3) * size)
    context.setVelocitiesToTemperature(temperature)

    mm.LocalEnergyMinimizer.minimize(context, 10, 100)

    # rudimentary T* basically
    # -1 -> not bonded
    # id of neighbor -> bonded to neighbor
    # -2 -> bonded to two neighbors
    active = np.repeat(-1, n_particles)

    # simulation
    sim_time = context.getTime()

    print("Start.", file=stderr)
    frame_count = 0

    start_real_time = time.time_ns()

    while sim_time < end_time:
        frame_count += 1
        state = context.getState(getPositions=True)
        sim_time = state.getTime()
        pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)

        # bonding algorithm
        if detection_modification(bonds, angles, active, pos):
            context.reinitialize(preserveState=True)
        # writing output
        print_frame(stdout, state, active)

        # n fully bonded
        n_free = np.count_nonzero(active == -1)
        n_full = np.count_nonzero(active == -2)
        n_edge = np.count_nonzero(active >= 0)
        perc = 100. * sim_time / end_time

        stderr.write(f"\rTime: {sim_time.value_in_unit(picosecond):.1f}" +
                     f" fs ({perc:.1f}%)\tParticles free: {n_free} " +
                     f"edge: {n_edge} full: {n_full}\033[K")
        integrator.step(steps_per_frame)

    end_real_time = time.time_ns()
    real_time_ms = (end_real_time - start_real_time) / 1000000
    print(f"\nDone. {frame_count} frames and {frame_count * steps_per_frame}" +
          f" steps in {real_time_ms:.0f} ms.", file=stderr)


def print_frame(file, state, active):
    positions = state.getPositions()
    natoms = len(positions)
    time = state.getTime() * 1000.0
    file.write(f"{natoms}\n")
    file.write(f"Bend or Break, t={time.value_in_unit(picosecond):.1f}\n")
    for i, pos in enumerate(positions):
        name = "N"
        if active[i] == -1:
            # free
            name = "C"
        elif active[i] == -2:
            # completely bonded
            name = "B"
        file.write(f"{name} {pos.x % size} {pos.y % size} {pos.z % size}\n")


if __name__ == "__main__":
    random.seed()
    simulate()
    
