from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stderr, stdout


def simulate():
    steps_per_frame = 100

    mass = 10.0
    charge = 1.0
    sigma = 0.1  # nm
    epsilon = 1.0  # kJ
    temperature = 300 * kelvin
    friction = 100 / picosecond
    step_size = 0.005 * picoseconds
    size = 10.0  # nm
    end_time = 100.0 * picoseconds

    system = System()

    system.addParticle(mass)
    system.addParticle(mass)

    nonbond = NonbondedForce()
    system.addForce(nonbond)
    nonbond.addParticle(charge, sigma, epsilon)
    nonbond.addParticle(charge, sigma, epsilon)

    bonds = HarmonicBondForce()
    bonds_idx = system.addForce(bonds)
    bonds.addBond(0, 1, 1.0, 100.0)  # FIXME: Hardcoding.

    integrator = LangevinMiddleIntegrator(temperature, friction, step_size)
    context = Context(system, integrator)

    half_size = size / 2
    initial_positions = [
        Vec3(half_size - sigma, half_size, half_size),
        Vec3(half_size + sigma, half_size, half_size),
    ]
    initial_velocities = [Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0)]
    context.setPositions(initial_positions)
    context.setVelocities(initial_velocities)

    context.setPeriodicBoxVectors(
        Vec3(size, 0.0, 0.0), Vec3(0.0, size, 0.0), Vec3(0.0, 0.0, size)
    )  # FIXME: Is there a nicer constructor here?

    # simulation.minimizeEnergy()

    time = context.getTime()

    while time < end_time / 2:
        state = context.getState(getPositions=True)

        print_frame(stdout, state)
        time = state.getTime()

        print(f"time: {time.value_in_unit(picoseconds):.1f} fs", file=stderr)
        integrator.step(steps_per_frame)

    # We remove the bond now.
    system.removeForce(bonds_idx)
    new_bonds = HarmonicBondForce()
    system.addForce(new_bonds)

    context.reinitialize(preserveState=True)
    print("Removed the bond.", file=stderr)

    while time < end_time:
        state = context.getState(getPositions=True)

        print_frame(stdout, state)
        time = state.getTime()

        print(f"time: {time.value_in_unit(picoseconds):.1f} fs", file=stderr)
        integrator.step(steps_per_frame)

    print("Done.", file=stderr)


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


if __name__ == "__main__":
    simulate()
