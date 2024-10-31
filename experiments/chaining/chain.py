# OpenMM example of a liquid of n particles, where harmonic bonds can be formed
# if two beads come within a cutoff distance

from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stderr, stdout
import random
import math
import numpy as np
import time

# input config
steps_per_frame = 250  # do frames (= bond calculations) every n steps
mass = 10.0
charge = 0.0
sigma = 0.1  # nm
epsilon = 1.  # kJ
temperature = 300 * kelvin
friction = 10 / picosecond
step_size = 0.005 * picoseconds
size = 15.  # nm, period box
end_time = 2000.0 * picoseconds
bond_length = 0.25 # nm
bond_force = 100.0
bond_cutoff = 1.0 # nm, below this threshold a bond is formed
n_particles = 150

bond_angle = math.pi # 180 degrees
angle_force = 100.0

# bond force and integrator choices
def gen_force():
    f = HarmonicBondForce()
    f.setUsesPeriodicBoundaryConditions(True)
    return f

def add_force(force, id0, id1):
    force.addBond(id0, id1, bond_length, bond_force)

def gen_angle_force():
    f = HarmonicAngleForce()
    f.setUsesPeriodicBoundaryConditions(True)
    return f

def add_angle(force, id0, id1, id2):
    force.addAngle(id1, id0, id2, bond_angle, angle_force)
    
def get_integrator():
    return LangevinMiddleIntegrator(temperature, friction, step_size)    
#    return VerletIntegrator(step_size)

def pdist(v1, v2, size):
    # periodic boundary adjusted difference between two Vec3
    # based on https://github.com/openmm/openmm/blob/master/platforms/reference/src/SimTKReference/ReferenceForce.cpp getDeltaRPeriodic

    diff = v1 - v2
    base = np.floor(diff / size + 0.5) * size
    return np.sqrt(np.dot(diff-base, diff-base))

assert(np.abs(pdist(np.array([0.,0.,1.]), np.array([0.,0.,4.]), 5.) - 2.) < 0.00001)
assert(np.abs(pdist(np.array([0.,0.,1.]), np.array([1.,0.,0.]), 5.) - np.sqrt(2)) < 0.00001)

def simulate():

    # context pre-setup, except particles
    system = System()
    nonbond = NonbondedForce()
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
    context = Context(system, integrator)
    context.setPeriodicBoxVectors(
        Vec3(size, 0.0, 0.0), Vec3(0.0, size, 0.0), Vec3(0.0, 0.0, size)
    )
    context.setPositions(np.random.rand(n_particles, 3) * size)
    context.setVelocitiesToTemperature(temperature)

    LocalEnergyMinimizer.minimize(context, 10, 100)
    
    # bonding setup
    # -1 -> not bonded
    # id of neighbor -> bonded to neighbor
    # -2 -> bonded to two neighbors
    active = np.repeat(-1, n_particles)
 
    def make_bond(i0, i1):
        add_force(bonds, i0, i1)
        if (active[i0] == -1):
            active[i0] = i1
        else:
            add_angle(angles, i0, active[i0], i1)
            active[i0] = -2
        if (active[i1] == -1):
            active[i1] = i0
        else:
            add_angle(angles, i1, active[i1], i0)
            active[i1] = -2

    # simulation
    sim_time = context.getTime()

    print("Start.", file=stderr)
    frame_count = 0

    start_real_time = time.time_ns()
    
    while sim_time < end_time:
        frame_count += 1
        state = context.getState(getPositions=True)
        sim_time = state.getTime()
        pos = state.getPositions(asNumpy=True).value_in_unit(nanometers)
        # bonding algorithm
        state_changes = False
        for i in range(len(active)):
            if active[i] == -2:
                continue
            neighbor = active[i]
            for j in range(i):
                if active[j] == -2:
                    continue
                new_neighbor = active[j]
                if i == j or neighbor == j or \
                        (neighbor >= 0 and neighbor == new_neighbor):
                    # no self reactions, no neighbor reactions, no neighbors neighbor reactions
                    continue
                
                dist = pdist(pos[i], pos[j], size)
                if dist < bond_cutoff:
                    make_bond(i, j)
                    state_changes = True
                    break
        
        if state_changes:
            context.reinitialize(preserveState=True)
            
        # writing output
        print_frame(stdout, state, active)
    
        # n fully bonded
        n_free = np.count_nonzero(active == -1)
        n_full = np.count_nonzero(active == -2)
        n_edge = np.count_nonzero(active >= 0)
        perc = 100. * sim_time / end_time

        stderr.write(f"\rTime: {sim_time.value_in_unit(picoseconds):.1f} fs ({perc:.1f}%)\tParticles free: {n_free} edge: {n_edge} full: {n_full}\033[K")
        integrator.step(steps_per_frame)


    end_real_time = time.time_ns()
    real_time_ms = (end_real_time - start_real_time) / 1000000
    print(f"\nDone. {frame_count} frames and {frame_count * steps_per_frame} steps in {real_time_ms:.0f} ms.", file=stderr)


def print_frame(file, state, active):
    positions = state.getPositions()
    natoms = len(positions)
    time = state.getTime() * 1000.0
    file.write(f"{natoms}\n")
    # Write the title. Reference a Polvo song.
    file.write(f"Bend or Break, t={time.value_in_unit(picoseconds):.1f}\n")
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
