#!/usr/bin/env python3

import openmm as mm
import openmm.app as mmapp
from openmm.unit import nanometer, picosecond, md_unit_system
from .utils import backup_try
from collections import OrderedDict
import numpy as np

from .forces.combined_bending_torsion import CombinedBendingTorsion
from .forces.force import Force
from .forces.g96angle import G96Angle
from .forces.harmonic_angle import HarmonicAngle
from .forces.harmonic_bond import HarmonicBond
from .forces.cubic_bond import CubicBond
from .forces.morse_bond import MorseBond
from .forces.fene_bond import FENEBond
from .forces.improper_dihedral import ImproperDihedral
from .forces.nonbonded import NonBonded, ExclusionHelper
from .forces.proper_dihedral import ProperDihedral
from .forces.rbtorsion import RBTorsion
from .forces.restricted_angle import RestrictedAngle
from .vsites.two_fd import VSite2fd
from .vsites.three_fd import VSite3fd
from .vsites.three_fad import VSite3fad
from .vsites.three_out import VSite3out
from .vsites.four_fdn import VSite4fdn
from .vsites.weighed_average import VSiteWeighedAverage
from .vsites.center_of_mass import VSiteCenterOfMass


class SysStar():

    _system: mm.System
    _context: mm.Context
    _reinitialize: bool
    _integrator: mm.Integrator
    nonbonded_force: NonBonded
    exclusions: ExclusionHelper
    _forces_list: list[mm.Force]  # to keep track of indices
    context_initialized: bool
    _part_list: list[(str, str, float, float)]

    """All the constraints in the system
    indices are called constraint_id
    values are (i: part_id, j: part_id, length: float)
    """
    _constraint_list: list[(int, int, float)]

    """Atom types to look up default charges and default masses
    values are atom_type: string, (mass: float, charge: float)
    """
    _atom_types: dict[str, (float, float)]

    modular_forces: list[Force]

    def __init__(self, epsilon_r, nonbonded_cutoff):
        self.epsilon_r = epsilon_r
        self.nonbonded_cutoff = nonbonded_cutoff
        self._part_list = []
        self._constraint_list = []
        self._atom_types = OrderedDict()
        self.context_initialized = False
        self._reinitialize = True
        self._forces_list = []
        self.harmonic_bond = HarmonicBond(self)
        self.morse_bond = MorseBond(self)
        self.cubic_bond = CubicBond(self)
        self.fene_bond = FENEBond(self)
        self.harmonic_angle = HarmonicAngle(self)
        self.proper_dihedral = ProperDihedral(self)
        self.improper_dihedral = ImproperDihedral(self)
        self.g96_angle = G96Angle(self)
        self.restricted_angle = RestrictedAngle(self)
        self.combined_bending_torsion = CombinedBendingTorsion(self)
        self.rb_torsion = RBTorsion(self)
        self.vsite_2fd = VSite2fd(self)
        self.vsite_3fd = VSite3fd(self)
        self.vsite_3fad = VSite3fad(self)
        self.vsite_3out = VSite3out(self)
        self.vsite_4fdn = VSite4fdn(self)
        self.vsite_avg = VSiteWeighedAverage(self)
        self.vsite_com = VSiteCenterOfMass(self)
        self.nonbonded_force = NonBonded(self)
        self.exclusions = self.nonbonded_force.get_exclusion_helper()
        self.modular_forces = [
            self.harmonic_bond, self.harmonic_angle,
            self.proper_dihedral, self.improper_dihedral,
            self.g96_angle, self.restricted_angle,
            self.combined_bending_torsion, self.rb_torsion,
            self.morse_bond, self.cubic_bond, self.fene_bond,
            self.vsite_2fd, self.vsite_3fd, self.vsite_4fdn,
            self.vsite_avg, self.vsite_3fad, self.vsite_3out, self.vsite_com,
            self.nonbonded_force, self.exclusions
        ]

        self.vsites = []

    def get_defaults(self, part_type, charge, mass):
        defaults = self._atom_types.get(part_type)
        if defaults is None:
            raise ValueError("Attempt to add particle with unknown atom type"
                             f"Particle type: {part_type}.")
        dcharge, dmass = defaults
        return (charge if charge is not None else dcharge,
                mass if mass is not None else dmass)

    def add_particle(self, part_name, part_type, charge, mass):
        """Adds a particle to the list, and returns its part_id"""
        if self.context_initialized:
            raise ValueError("Do this operation before context is initialized")
        charge, mass = self.get_defaults(part_type, charge, mass)
        self._part_list.append((part_name, part_type, charge, mass))
        return len(self._part_list) - 1

    def update_particle(self, part_id, part_name, part_type, charge, mass):
        old_name, old_type, old_charge, old_mass = self._part_list[part_id]
        charge, mass = self.get_defaults(part_type, charge, mass)
        self._part_list[part_id] = (part_name, part_type, charge, mass)
        if not self.context_initialized:
            return
        if old_type != part_type or old_charge != charge:
            # Change of LJ params plus self correction force (called by nb)
            self.nonbonded_force.update_params(
                part_id, part_type, charge, old_charge != charge
            )

        if old_mass != mass:
            self._system.setParticleMass(part_id, mass)
        self._reinitialize = True

    def build_system(self):
        self._system = mm.System()
        for (_, _, _, mass) in filter(None, self._part_list):
            self._system.addParticle(mass)
        for (i, j, length) in filter(None, self._constraint_list):
            self._system.addConstraint(i, j, length)

    def get_particle_name_type(self, i):
        name, type, _, _ = self._part_list[i]
        return name, type

    def get_particle_details(self, i):
        """Returns the particle's name, type, charge, mass"""
        return self._part_list[i]

    def len_particles(self):
        return len(self._part_list)

    def add_constraint(self, i, j, length):
        """Adds a constraint to the list, returns its constraint_id"""
        self._constraint_list.append((i, j, length))
        if self.context_initialized:
            raise Exception("Can't add constraint during run")
        return len(self._constraint_list) - 1

    def get_constraint_members(self, constraint_id):
        i, j, _ = self._constraint_list[constraint_id]
        return (i, j)

    def remove_constraint(self, constraint_id):
        raise NotImplementedError("Constraints cannot be safely removed")

    def add_atom_type(self, type, charge, mass):
        if self.context_initialized:
            raise ValueError("Cannot do this after context is initialized")
        self._atom_types[type] = (charge, mass)

    def add_nb_type(self, type1, type2, V, W):
        if self.context_initialized:
            raise ValueError("Cannot do this after context is initialized")
        self.nonbonded_force._nb_types[(type1, type2)] = (V, W)

    def build_context(self, integrator, periodicBoxVectors, platform=None):
        """context_initialized flips the state of S* in a way
        the parsing of a topology and the addition of all particles, bonds, ...
        should happen before calling build_context
        """
        if self.context_initialized:
            raise ValueError("Cannot do this after context is initialized")
        # Build system+forces
        self.build_system()
        for force in self._forces_list:
            self._system.addForce(force)
        self.context_initialized = True
        self._reinitialize = False
        self._integrator = integrator
        for modular_force in self.modular_forces:
            modular_force.build()
        # Build context
        self._periodic_box = periodicBoxVectors
        self._system.setDefaultPeriodicBoxVectors(*periodicBoxVectors)
        if platform is None:
            self._context = mm.Context(self._system, integrator)
        else:
            self._context = mm.Context(self._system, integrator, platform)
        self._context.setPeriodicBoxVectors(*periodicBoxVectors)

    def reinitialize(self):
        if not self.context_initialized:
            raise Exception("Initialize the context first")
        for modular_force in self.modular_forces:
            modular_force.build()
        if self._reinitialize:
            self._context.reinitialize(preserveState=True)

    def set_positions(self, positions):
        if not self.context_initialized:
            raise Exception("Initialize the context first")
        self._context.setPositions(positions)

    def generate_velocities(self, temp):
        if not self.context_initialized:
            raise Exception("Initialize the context first")
        self._context.setVelocitiesToTemperature(temp)

    def minimize_energy(self, tolerance=10, max_steps=0):
        # 0 max_steps => until converged
        if not self.context_initialized:
            raise Exception("Initialize the context first")
        mm.LocalEnergyMinimizer.minimize(self._context, tolerance, max_steps)

    def add_force(self, force):
        """Adds a force to S*
        returns if caller should call reinitialize"""
        if self.context_initialized:
            self._forces_list.append(force)
            self._system.addForce(force)
            self._reinitialize = True
            return True
        else:
            # contents of _forces_list automatically get added to the system
            # during self.build_context()
            self._forces_list.append(force)
            return False

    def do_steps(self, steps):
        if not self.context_initialized:
            raise Exception("Initialize the context first")
        self._integrator.step(steps)
        if self._xtc is not None:
            self._xtc.interval = steps
            pos, _ = self.get_positions()
            self._xtc.writeModel(pos)

    def get_positions(self):
        if not self.context_initialized:
            raise Exception("Initialize the context first")
        # get the state without enforcePeriodicBox
        # since enforcePeriodicBox really doesn't like bonds formed across
        # boundaries
        state = self._context.getState(positions=True)
        box = state.getPeriodicBoxVectors()
        box_x = box[0].x
        box_y = box[1].y
        box_z = box[2].z
        # TODO: don't assume 90 degree angles in martini_daemon in general
        # these asserts are there to ensure only 90 degree boxes are ran
        assert box[1].x == 0.
        assert box[2].x == 0.
        assert box[2].y == 0.
        pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)
        pos[:, 0] = np.remainder(pos[:, 0], box_x)
        pos[:, 1] = np.remainder(pos[:, 1], box_y)
        pos[:, 2] = np.remainder(pos[:, 2], box_z)
        return pos, np.array([box_x, box_y, box_z])

    def get_state(self):
        """Returns a state with forces and energies, used in test.py."""
        if not self.context_initialized:
            raise Exception("Initialize the context first")
        state = self._context.getState(positions=True, velocities=True,
                                       forces=True, energy=True)
        return state

    def apply_constraints(self):
        if not self.context_initialized:
            raise Exception("Initialize the context first")
        self._context.applyConstraints(tol=1e-10)

    _xtc: mmapp.XTCFile = None

    def set_xtc_path(self, path):
        if not self.context_initialized:
            raise Exception("Initialize the context first")
        backup_try(path)
        mmtopol = mmapp.Topology()
        mmtopol._numAtoms = self.len_particles()
        mmtopol._periodicBoxVectors = self._periodic_box
        timestep = self._integrator.getStepSize()
        self._xtc = mmapp.XTCFile(path, mmtopol, timestep)

    def write_gro(self, path):
        if not self.context_initialized:
            raise Exception("Initialize the context first")
        backup_try(path)
        state = self._context.getState(positions=True, velocities=True)
        pos, _ = self.get_positions()
        vel = state.getVelocities(asNumpy=True).\
            value_in_unit_system(md_unit_system)
        natoms = len(pos)
        time = state.getTime()
        with open(path, "w") as file:
            file.write(f"t={time.value_in_unit(picosecond):.1f} ps\n")
            file.write(f"{natoms}\n")
            for i in range(natoms):
                cpos = pos[i]
                cvel = vel[i]
                atom_name = self._part_list[i][0]
                atom_index = i + 1
                # TODO resname, resid
                file.write(f"{atom_index:5}{atom_name:5}{atom_name:>5}"
                           f"{atom_index:5}{cpos[0]:8.3f}{cpos[1]:8.3f}"
                           f"{cpos[2]:8.3f}{cvel[0]:8.4f}{cvel[1]:8.4f}"
                           f"{cvel[2]:8.4f}\n")
            v1, v2, v3 = (v.value_in_unit(nanometer)
                          for v in state.getPeriodicBoxVectors())
            file.write(
                f"{v1[0]:.4f} {v2[1]:.4f} {v3[2]:.4f} {v1[1]:.4f} {v1[2]:.4f} "
                f"{v2[0]:.4f} {v2[2]:.4f} {v3[0]:.4f} {v3[1]:.4f}\n"
            )

    def dump(self):
        backup_try("sys.dump")
        with open("sys.dump", "w") as file:
            print("==== SysStar Dump ====", file=file)
            print("Particles:", self._part_list, file=file)
            for force in self.modular_forces:
                print(force.__class__.__name__, force._list, file=file)
            print("Exclusions:", self._exclusion_list, file=file)
            print("Constraints:", self._constraint_list, file=file)
