from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stderr, stdout
import martini_openmm as martini

def print_frame(file, state):
    positions = state.getPositions()

    natoms = len(positions)
    time = state.getTime() * 1000.0  # fs to ps.

    # Write the title. Reference a Polvo song.
    print(f"Bend or Break, t={time.value_in_unit(picoseconds):.1f}", file=file)
    # Write the number of atoms in the system.
    print(f"{natoms}", file=file)

    for i, pos in enumerate(positions):
        # Write the atom line.
        print(
            f"{i:5}{'DUMMY':>5}{'A':5}{i:5}{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}",
            file=file,
        )

    # Write the box vectors.
    v1, v2, v3 = (v.value_in_unit(nanometers) for v in state.getPeriodicBoxVectors())
    # v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
    print(
        f"{v1[0]:.4f} {v2[1]:.4f} {v3[2]:.4f} {v1[1]:.4f} {v1[2]:.4f} {v2[0]:.4f} {v2[2]:.4f} {v3[0]:.4f} {v3[1]:.4f}",
        file=file,
    )

end_time = 100.0 * picoseconds
steps_per_frame = 100

# Generating System using Martini OpenMM
conf = GromacsGroFile("system.gro")
box_vectors = conf.getPeriodicBoxVectors()

top  = martini.MartiniTopFile("topol.top",
                periodicBoxVectors=box_vectors,
                defines={},
                epsilon_r=15)
system = top.create_system(nonbonded_cutoff=1.1 * nanometer)
# Continuing as normal with system setup
integrator = LangevinIntegrator(310 * kelvin, 10.0 / picosecond, 20 * femtosecond)
integrator.setRandomNumberSeed(0)
context  = Context(system, integrator)
context.setPositions(conf.getPositions())

time = context.getTime()
while time < end_time / 4:
    state = context.getState(getPositions=True)
    print_frame(stdout, state)
    time = state.getTime()
    print(f"time: {time.value_in_unit(picoseconds):.1f} fs", file=stderr)
    integrator.step(steps_per_frame)

# Getting information about the bond
force = system.getForce(4)
bond_forces = [force.getBondParameters(idx) for idx in range(force.getNumBonds())]
# Removing bond
system.removeForce(4)
bond_idx = system.addForce(HarmonicBondForce())

context.reinitialize(preserveState=True)
print("Removed the bond.", file=stderr)
while time < end_time / 2:
    state = context.getState(getPositions=True)
    print_frame(stdout, state)
    time = state.getTime()
    print(f"time: {time.value_in_unit(picoseconds):.1f} fs", file=stderr)
    integrator.step(steps_per_frame)

# Adding bond back in
force = system.getForce(bond_idx)
for p1, p2, *a in bond_forces:
    force.addBond(p1, p2, *a)

context.reinitialize(preserveState=True)
print("Add bond.", file=stderr)
while time < end_time * 0.75:
    state = context.getState(getPositions=True)
    print_frame(stdout, state)
    time = state.getTime()
    print(f"time: {time.value_in_unit(picoseconds):.1f} fs", file=stderr)
    integrator.step(steps_per_frame)

# Removing bond again
system.removeForce(bond_idx)

context.reinitialize(preserveState=True)
print("Removed the bond.", file=stderr)
while time < end_time :
    state = context.getState(getPositions=True)
    print_frame(stdout, state)
    time = state.getTime()
    print(f"time: {time.value_in_unit(picoseconds):.1f} fs", file=stderr)
    integrator.step(steps_per_frame)


print("Done.", file=stderr)


