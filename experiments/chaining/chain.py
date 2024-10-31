# OpenMM example of a liquid of n particles, where harmonic bonds can be formed
# if two beads come within a cutoff distance

from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stderr, stdout
import random
import math

# input config
steps_per_frame = 100
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
angle_force = 10.0

# bond force and integrator choices
def gen_force():
    return HarmonicBondForce()

def add_force(force, id0, id1):
    force.addBond(id0, id1, bond_length, bond_force)

def gen_angle_force():
    return HarmonicAngleForce()

def add_angle(force, id0, id1, id2):
    force.addAngle(id1, id0, id2, bond_angle, angle_force)
    
def get_integrator():
    return LangevinMiddleIntegrator(temperature, friction, step_size)    
#    return VerletIntegrator(step_size)

def rand_vec3(min, max):
    # generates a random Vec3 within the sim box

    return Vec3(random.uniform(min, max),
                random.uniform(min, max),
                random.uniform(min, max))

def pdiff(x1, x2, size):
    # periodic boundary adjusted difference between two floats
    # based on https://github.com/openmm/openmm/blob/master/platforms/reference/src/SimTKReference/ReferenceForce.cpp periodicDifference

    diff = x1 - x2
    base = math.floor(diff / size + 0.5) * size
    return diff - base

def pdist(v1, v2, size):
    # periodic boundary adjusted difference between two Vec3
    # based on https://github.com/openmm/openmm/blob/master/platforms/reference/src/SimTKReference/ReferenceForce.cpp getDeltaRPeriodic

    diff = Vec3(
        pdiff(v1.x, v2.x, size),
        pdiff(v1.y, v2.y, size),
        pdiff(v1.z, v2.z, size)
    )

    return math.sqrt(diff.x ** 2 + diff.y ** 2 + diff.z ** 2)

assert(pdiff(1., 2., 5.) - 1. < 0.00001)
assert(pdiff(1., 4., 5.) - 2. < 0.00001)
assert(pdist(Vec3(0.,0.,1.), Vec3(1.,0.,0.), 5.) - math.sqrt(2) < 0.00001)

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
    initial_positions = []  # list of Vec3, nm
    for i in range(n_particles):
        system.addParticle(mass)
        nonbond.addParticle(charge, sigma, epsilon)
        initial_positions.append(rand_vec3(0., size))
    

    # context setup
    context = Context(system, integrator)
    context.setPeriodicBoxVectors(
        Vec3(size, 0.0, 0.0), Vec3(0.0, size, 0.0), Vec3(0.0, 0.0, size)
    )
    context.setPositions(initial_positions)
    context.setVelocitiesToTemperature(temperature)

    LocalEnergyMinimizer.minimize(context, 10, 100)
    
    # bonding setup
    # list of reactive sites
    # of the form (id, neighbor_id)
    # if neighbor_id is 0 then it's not bonded yet
    # if two bonds are fulfilled, it is removed from the reactive list by setting the id to -1
    reactive = [] # TODO: linked list is a better data structure for this
    n_reactive = 0
    for i in range(n_particles):
        reactive.append([i, -1])
        n_reactive += 1

    def make_bond(i0, i1):
        add_force(bonds, i0, i1)

    # simulation
    time = context.getTime()

    print("Start.", file=stderr)
    while time < end_time:
        state = context.getState(getPositions=True)
        pos = state.getPositions()
        # bonding algorithm
        state_changes = False
        for i in range(n_reactive):
            id_this = reactive[i][0]
            if id_this == -1:
                continue
            id_neighbor = reactive[i][1]
            for j in range(n_reactive):
                id_new = reactive[j][0]
                if id_new == -1:
                    continue
                id_new_neighbor = reactive[j][1]
                if id_this == id_new or id_neighbor == id_new or \
                        (id_neighbor >= 0 and id_neighbor == id_new_neighbor):
                    # no self reactions, no neighbor reactions, no neighbors neighbor reactions
                    continue
                
                dist = pdist(pos[id_this], pos[id_new], size)
                if dist < bond_cutoff:
                    make_bond(id_this, id_new)
                    if id_neighbor == -1:
                        reactive[id_this][1] = id_new
                    else:
                        reactive[id_this][0] = -1
                        add_angle(angles, id_this, id_neighbor, id_new)
                        n_reactive -= 1
                    if id_new_neighbor == -1:
                        reactive[id_new][1] = id_this
                    else:
                        reactive[id_new][0] = -1
                        add_angle(angles, id_new, id_this, id_new_neighbor)
                        n_reactive -= 1
                    state_changes = True
                    break
        

        if state_changes:
            context.reinitialize(preserveState=True)
            
                

        # writing output
        print_frame(stdout, state)
        time = state.getTime()
        
        # n fully bonded
        n_bonded = n_particles - n_reactive
        perc = 100. * time / end_time
        perc_bonded = 100. * n_bonded / n_particles
        stderr.write(f"\rTime: {time.value_in_unit(picoseconds):.1f} fs ({perc:.1f}%)\t\tBonded: {n_bonded} ({perc_bonded:.1f}%)")
        integrator.step(steps_per_frame)


    print("\nDone.", file=stderr)


def print_frame(file, state):
    positions = state.getPositions()

    natoms = len(positions)
    time = state.getTime() * 1000.0  # fs to ps.

    # Write the title. Reference a Polvo song.
    print(f"Bend or Break, t={time.value_in_unit(picoseconds):.1f}", file=file)
    # Write the number of atoms in the system.
    print(f"{natoms}", file=file)

    for i, pos in enumerate(positions):
        # Write the atom line. Modulo with box size to display position within pbc
        name = "DUMMY"
        elem = "A"
        velx = 0.0
        print(
            f"{i:5}{name:>5}{elem:5}{i:5}{(pos.x % size):8.3f}{(pos.y % size):8.3f}{(pos.z % size):8.3f}{velx:8.4f}{0.:8.4f}{0.:8.4f}",
            file=file,
        )

    # Write the box vectors.
    v1, v2, v3 = (v.value_in_unit(nanometers) for v in state.getPeriodicBoxVectors())
    # v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
    print(
        f"{v1[0]:.4f} {v2[1]:.4f} {v3[2]:.4f} {v1[1]:.4f} {v1[2]:.4f} {v2[0]:.4f} {v2[2]:.4f} {v3[0]:.4f} {v3[1]:.4f}",
        file=file,
    )


if __name__ == "__main__":
    random.seed()
    simulate()
