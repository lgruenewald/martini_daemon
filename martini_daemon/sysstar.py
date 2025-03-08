#!/usr/bin/env python3

import openmm as mm
import openmm.app as mmapp
from openmm.unit import nanometer, picosecond, md_unit_system
from .utils import backup_try
from collections import OrderedDict
import numpy as np
import sys

from .forces.combined_bending_torsion import CombinedBendingTorsion
from .forces.constraint import Constraint
from .forces.force import Force
from .forces.g96angle import G96Angle
from .forces.harmonic_angle import HarmonicAngle
from .forces.harmonic_bond import HarmonicBond
from .forces.cubic_bond import CubicBond
from .forces.morse_bond import MorseBond
from .forces.fene_bond import FENEBond
from .forces.improper_dihedral import ImproperDihedral
from .forces.nonbonded import NonBonded, ExclusionHelper
from .forces.pairs import Pairs
from .forces.proper_dihedral import ProperDihedral
from .forces.rbtorsion import RBTorsion
from .forces.restricted_angle import RestrictedAngle
from .forces.restricted_dihedral import RestrictedDihedral
from .forces.g96bond import G96Bond
from .forces.connection import Connection
from .forces.distance_restraint import DistanceRestraint
from .forces.urey_bradley import UreyBradley
from .forces.cross_bond_bond import CrossBondBond
from .forces.cross_bond_angle import CrossBondAngle
from .forces.quartic_angle import QuarticAngle
from .forces.linear_angle import LinearAngle
from .vsites.two_fd import VSite2fd
from .vsites.three_fd import VSite3fd
from .vsites.three_fad import VSite3fad
from .vsites.three_out import VSite3out
from .vsites.four_fdn import VSite4fdn
from .vsites.weighed_average import VSiteWeighedAverage
from .vsites.center_of_mass import VSiteCenterOfMass
from .reporters.reporter import Reporter


class SysStar():

    _system: mm.System
    _context: mm.Context
    _reinitialize: bool
    _integrator: mm.Integrator
    nonbonded_force: NonBonded
    exclusions: ExclusionHelper
    _forces_list: list[mm.Force]  # to keep track of indices
    context_initialized: bool
    # name, resid, resname, type, charge, mass
    _part_list: list[(str, int, str, str, float, float)]

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

    """Extra reporters that can write stuff to files every step / at the end
    """
    reporters: list[Reporter]

    def __init__(self, logger, epsilon_r, nonbonded_cutoff):
        self.epsilon_r = epsilon_r
        self.nonbonded_cutoff = nonbonded_cutoff
        self._part_list = []
        self._constraint_list = []
        self._atom_types = OrderedDict()
        self.context_initialized = False
        self._reinitialize = True
        self._forces_list = []
        self.constraint = Constraint(self)
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
        self.g96_bond = G96Bond(self)
        self.connection = Connection(self)
        self.distance_restraint = DistanceRestraint(self)
        self.urey_bradley = UreyBradley(self)
        self.restricted_dihedral = RestrictedDihedral(self)
        self.cross_bond_bond = CrossBondBond(self)
        self.cross_bond_angle = CrossBondAngle(self)
        self.quartic_angle = QuarticAngle(self)
        self.linear_angle = LinearAngle(self)
        self.pairs = Pairs(self)
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
            self.constraint,
            self.harmonic_bond, self.harmonic_angle,
            self.proper_dihedral, self.improper_dihedral,
            self.g96_angle, self.restricted_angle,
            self.combined_bending_torsion, self.rb_torsion,
            self.morse_bond, self.cubic_bond, self.fene_bond,
            self.g96_bond, self.connection, self.distance_restraint,
            self.urey_bradley, self.restricted_dihedral,
            self.cross_bond_bond, self.cross_bond_angle,
            self.quartic_angle, self.linear_angle,
            self.vsite_2fd, self.vsite_3fd, self.vsite_4fdn,
            self.vsite_avg, self.vsite_3fad, self.vsite_3out, self.vsite_com,
            self.nonbonded_force, self.exclusions, self.pairs
        ]

        self.vsites = []

        self.reporters = []

        self.last_resid = 0
        self.logger = logger

    def get_defaults(self, part_type, charge, mass):
        defaults = self._atom_types.get(part_type)
        if defaults is None:
            raise ValueError("Attempt to add particle with unknown atom type"
                             f"Particle type: {part_type}.")
        dcharge, dmass = defaults
        return (charge if charge is not None else dcharge,
                mass if mass is not None else dmass)

    def new_residue(self):
        self.last_resid += 1

    def add_particle(self, part_name, resname, part_type, charge, mass):
        """Adds a particle to the list, and returns its part_id"""
        if self.context_initialized:
            raise ValueError("Do this operation before context is initialized")
        charge, mass = self.get_defaults(part_type, charge, mass)
        self._part_list.append((part_name, self.last_resid, resname, part_type, charge, mass))
        return len(self._part_list) - 1

    def rename(self, part_id, new_name):
        self._part_list[part_id][0] = new_name

    def retype(self, part_id, new_type):
        _, _, _, old_type, charge, _ = self._part_list[part_id]
        self._part_list[part_id][3] = new_type
        if self.context_initialized and old_type != new_type:
            self.nonbonded_force.update_params(
                part_id, new_type, charge, False
            )

    def recharge(self, part_id, new_charge):
        _, _, _, type, old_charge, _ = self._part_list[part_id]
        self._part_list[part_id][4] = new_charge
        if self.context_initialized and old_charge != new_charge:
            self.nonbonded_force.update_params(
                part_id, type, new_charge, True
            )

    def remass(self, part_id, new_mass):
        old_mass = self._part_list[part_id][5]
        self._part_list[part_id][5] = new_mass
        if self.context_initialized and new_mass != old_mass:
            self._system.setParticleMass(part_id, new_mass)
            self._reinitialize = True

    def update_particle(self, part_id, part_type, charge, mass):
        self.retype(part_id, part_type)
        self.recharge(part_id, charge)
        self.remass(part_id, mass)

    def build_system(self):
        self._system = mm.System()
        for (_, _, _, _, _, mass) in filter(None, self._part_list):
            self._system.addParticle(mass)

    def get_particle_name_type(self, i):
        name, _, _, type, _, _ = self._part_list[i]
        return name, type

    def get_particle_details(self, i):
        """Returns the particle's type, charge, mass"""
        _, _, _, type, charge, mass = self._part_list[i]
        return (type, charge, mass)

    def len_particles(self):
        return len(self._part_list)

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
        self._integrator = integrator
        for modular_force in self.modular_forces:
            modular_force.build()
        self._reinitialize = False
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
        if len(positions) != len(self._part_list):
            raise ValueError(
                "Wrong length of coordinate file."
                f"Supplied {len(positions)} coordinates."
                f"Particle list contains {len(self._part_list)} particles."
                "Please check your .gro file."
            )
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

    def add_reporter(self, reporter):
        self.reporters.append(reporter)

    def write_xtc_frame(self, interval):
        pos, box = self.get_positions()
        self._xtc.interval = interval
        self._xtc.writeModel(pos)
        for reporter in self.reporters:
            reporter.on_xtc_frame(pos, box, self._xtc_name)

    def do_steps(self, steps):
        if not self.context_initialized:
            raise Exception("Initialize the context first")
        self._integrator.step(steps)

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
        # only works for 90 degree pbc
        assert box[1].x == 0.
        assert box[2].x == 0.
        assert box[2].y == 0.
        pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)
        if np.any(np.abs(pos) > 2147483.0):
            # would be too large to store without remaindering, so likely
            # the system blew up
            self.logger.error(
                "Particle coordinates too large, your system likely blew up. "
                f"Largest coordinate (abs value) is {np.max(np.abs(pos))}."
            )
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
        if len(path) < 5 or path[-4:] != ".xtc":
            raise Exception("Xtc path must end with .xtc")
        backup_try(path)
        for reporter in self.reporters:
            reporter.on_set_xtc_path(path[:-4])
        mmtopol = mmapp.Topology()
        mmtopol._numAtoms = self.len_particles()
        mmtopol._periodicBoxVectors = self._periodic_box
        timestep = self._integrator.getStepSize()
        self._xtc = mmapp.XTCFile(path, mmtopol, timestep)
        self._xtc_name = path[:-4]

    def write_gro(self, path):
        if not self.context_initialized:
            raise Exception("Initialize the context first")
        if len(path) < 5 or path[-4:] != ".gro":
            raise Exception("Gro path must end with .gro")
        backup_try(path)
        state = self._context.getState(positions=True, velocities=True)
        pos, box = self.get_positions()
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
                atom_name, resid, resname, *_ = self._part_list[i]
                atom_index = i + 1
                file.write(f"{resid:5}{resname:5}{atom_name:>5}"
                           f"{atom_index:5}{cpos[0]:8.3f}{cpos[1]:8.3f}"
                           f"{cpos[2]:8.3f}{cvel[0]:8.4f}{cvel[1]:8.4f}"
                           f"{cvel[2]:8.4f}\n")
            v1, v2, v3 = (v.value_in_unit(nanometer)
                          for v in state.getPeriodicBoxVectors())
            file.write(
                f"{v1[0]:.4f} {v2[1]:.4f} {v3[2]:.4f} {v1[1]:.4f} {v1[2]:.4f} "
                f"{v2[0]:.4f} {v2[2]:.4f} {v3[0]:.4f} {v3[1]:.4f}\n"
            )
        for reporter in self.reporters:
            reporter.on_write_gro(pos, box, path[:-4])

    def dump(self):
        backup_try("sys.dump")
        with open("sys.dump", "w") as file:
            print("==== SysStar Dump ====", file=file)
            print("Particles:", self._part_list, file=file)
            for force in self.modular_forces:
                print(force.__class__.__name__, force._list, file=file)
            print("Exclusions:", self._exclusion_list, file=file)
            print("Constraints:", self._constraint_list, file=file)
