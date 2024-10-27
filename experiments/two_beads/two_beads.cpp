/*
 * A small proof of concept two breaking the bond potential between two
 * particles in the middle of a simulation.
 *
 * Marieke Westendorp, 2024.
 */

#include "OpenMM.h"
#include "openmm/HarmonicBondForce.h"
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace OpenMM;

// Print a gro frame based on a `state`.
//
// https://manual.gromacs.org/archive/5.0.3/online/gro.html
void print_frame(std::FILE *stream, const State &state) {
  const auto positions = state.getPositions();

  const auto natoms = positions.size();
  const auto time = state.getTime() * 1000.0; // fs to ps.

  // Write the title. Reference a Polvo song.
  std::fprintf(stream, "Bend or Break, t=%.1f\n", time);
  // Write the number of atoms in the system.
  std::fprintf(stream, "%zu\n", natoms);

  for (std::vector<Vec3>::size_type i = 0; i < natoms; i++) {
    const auto pos = positions[i];

    // Write the atom line.
    std::fprintf(stream, "%5zu%-5s%5s%5zu%8.3f%8.3f%8.3f\n", i, "DUMMY", "A", i,
                 pos[0], pos[1], pos[2]);
  }

  // Write the box vectors.
  Vec3 v1, v2, v3 = {};
  state.getPeriodicBoxVectors(v1, v2, v3);
  // v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
  std::fprintf(stream, "%f %f %f %f %f %f %f %f %f\n", v1[0], v2[1], v3[2],
               v1[1], v1[2], v2[0], v2[2], v3[0], v3[1]);
}

void simulate() {
  const double mass = 10.0;
  const double charge = 1.0;
  const double sigma = 0.1;       //  nm
  const double epsilon = 1.0;     // kJ
  const double temperature = 300; // K
  const double friction = 100;    // per ps
  const double step_size = 0.005; // ps
  const double size = 10.0;       // nm

  System system;

  // Add two particles.
  system.addParticle(mass);
  system.addParticle(mass);

  // Add a nonbonded potential between the particles.
  NonbondedForce *nonbond = new NonbondedForce();
  system.addForce(nonbond);
  nonbond->addParticle(charge, sigma, epsilon);
  nonbond->addParticle(charge, sigma, epsilon);

  // Set up their bond.
  HarmonicBondForce *bonds = new HarmonicBondForce();
  const auto bonds_idx = system.addForce(bonds);
  bonds->addBond(0, 1, 1.0, 100.0);

  // No angles for now, since we have only two particles.
  // HarmonicAngleForce *angles = new HarmonicAngleForce();
  // system.addForce(angles);

  // Create an integrator and set up the context.
  LangevinMiddleIntegrator integrator(temperature, friction, step_size);
  Context context(system, integrator);

  // Set their initial positions.
  const double half_size = size / 2.0;
  std::vector<Vec3> initial_positions = {
      Vec3(half_size - sigma, half_size, half_size),
      Vec3(half_size + sigma, half_size, half_size)};
  std::vector<Vec3> initial_velocities = {Vec3(0.0, 0.0, 0.0),
                                          Vec3(0.0, 0.0, 0.0)};
  context.setPositions(initial_positions);
  context.setVelocities(initial_velocities);

  // Set the box vectors
  context.setPeriodicBoxVectors(Vec3(size, 0.0, 0.0), Vec3(0.0, size, 0.0),
                                Vec3(0.0, 0.0, size));

  // Time to simulate!
  const double end_time = 100.0;
  double time = context.getTime();
  const auto types = State::Positions;

  // First with the bonded potential in place.
  while (time < end_time / 2) {
    State state = context.getState(types);

    print_frame(stdout, state);
    time = state.getTime();

    std::fprintf(stderr, "time: %.1f fs\n", time);
    integrator.step(100);
  }

  // We remove the bond now.
  // This is where the interesting thing happens (in an admittedly boring way).
  system.removeForce(bonds_idx); // Remove the old bonds force.
  // Setting up a new bonds force description, but we don't add any bonds.
  //
  // Note that this is kind of redundant in our two particle system, because by
  // breaking the single bond that exists in the system, we may as well not
  // introduce the bond forces at all. But for larger systems where many obnds
  // will remain and only a single or very few bonds are broken, we would need
  // to set up the alternate bonds configuration like this.
  HarmonicBondForce *new_bonds = new HarmonicBondForce();
  system.addForce(new_bonds);
  // Don't add any bonds to the new bonds definition.

  // From the above, it is obvious to me that a more serious example or eventual
  // implementation of the principle should involve a custom object that keeps
  // track of the bonds and has associated functions for converting to the
  // OpenMM bonds force descriptions based on that bookkeeping.
  // In other words, this is just a silly little hack until serious bookkeeping
  // is introduced.

  context.reinitialize(true);
  std::fprintf(stderr, "Removed the bond.\n");

  // And we continue simulating the system with the bond removed.
  while (time < end_time) {
    State state = context.getState(types);

    print_frame(stdout, state);
    time = state.getTime();

    std::fprintf(stderr, "time: %.1f fs\n", time);
    integrator.step(100);
  }

  std::fprintf(stderr, "Done.\n");
}

int main() {
  try {
    simulate();
  } catch (const std::exception &e) {
    std::printf("ERROR: exception: %s\n", e.what());
    return 1;
  }
  return 0;
}
