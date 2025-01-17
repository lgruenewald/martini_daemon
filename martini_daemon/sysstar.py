#!/usr/bin/env python3

from __future__ import annotations
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
from .forces.improper_dihedral import ImproperDihedral
from .forces.proper_dihedral import ProperDihedral
from .forces.rbtorsion import RBTorsion
from .forces.restricted_angle import RestrictedAngle
from .vsites.three_fad import VSite3fad
from .vsites.three_out import VSite3out
from .vsites.weighed_average import VSiteWeighedAverage


class SysStar():

    _system: mm.System
    _context: mm.Context
    _reinitialize: bool
    _integrator: mm.Integrator
    _nb_force: mm.Force
    _nb_force_rebuild = False
    _es_self_correction_force: mm.Force
    _es_self_correction_force_rebuild = False
    _forces_list: list[mm.Force]  # to keep track of indices
    context_initialized: bool
    _part_list: list[(str, str, float, float)]
    _used_atom_types: OrderedDict[str, int]
    """All the exclusions in the system
    indices are caled excl_id
    values are (i: part_id, j: part_id)
    """
    _exclusion_list: list[(int, int)]
    """All the constraints in the system
    indices are called constraint_id
    values are (i: part_id, j: part_id, length: float)
    """
    _constraint_list: list[(int, int, float)]
    """Atom types to look up default charges and default masses
    values are atom_type: string, (mass: float, charge: float)
    """
    _atom_types: dict[str, (float, float)]
    """Non bonded parameters are for building the C6/C12 table
    values are (type1: string, type2: string), (V: float, W: float)
    """
    _nb_types: dict[(str, str), (float, float)]

    modular_forces: list[Force]

    def __init__(self, epsilon_r, nonbonded_cutoff):
        self._part_list = []
        self._used_atom_types = {}
        self._exclusion_list = []
        self._constraint_list = []
        self._atom_types = OrderedDict()
        self.context_initialized = False
        self._reinitialize = True
        self._nb_types = {}
        self._forces_list = []
        self.epsilon_r = epsilon_r
        self.nonbonded_cutoff = nonbonded_cutoff
        self.harmonic_bond = HarmonicBond(self)
        self.harmonic_angle = HarmonicAngle(self)
        self.proper_dihedral = ProperDihedral(self)
        self.improper_dihedral = ImproperDihedral(self)
        self.g96_angle = G96Angle(self)
        self.restricted_angle = RestrictedAngle(self)
        self.combined_bending_torsion = CombinedBendingTorsion(self)
        self.rb_torsion = RBTorsion(self)
        self.vsite_3fad = VSite3fad(self)
        self.vsite_3out = VSite3out(self)
        self.vsite_avg = VSiteWeighedAverage(self)
        self.modular_forces = [
            self.harmonic_bond, self.harmonic_angle,
            self.proper_dihedral, self.improper_dihedral,
            self.g96_angle, self.restricted_angle,
            self.combined_bending_torsion, self.rb_torsion,
            self.vsite_avg, self.vsite_3fad, self.vsite_3out
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
            # Change of LJ params plus self correction force
            part_type_id = self.use_atom_type(part_type)
            self._nb_force.setParticleParameters(part_id, [part_type_id, charge])
            if old_charge != charge:
                # TODO optimize this later maybe?
                self._es_self_correction_force_rebuild = True

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

    def use_atom_type(self, part_type):
        """Add new atom types for LJ, used in build_nb_force.

        It is not a problem to call this with atom types already present.

        Returns the atom type index
        """
        id = self._used_atom_types.get(part_type)
        if id is None:
            id = len(self._used_atom_types)
            self._used_atom_types[part_type] = id
            self._nb_force_rebuild = True
        return id

    def len_particles(self):
        return len(self._part_list)

    def build_nb_force(self):
        self._nb_force = mm.CustomNonbondedForce(
            "step(rcut-r)*(LJ - corr + ES);"
            "LJ = (W(type1, type2) / r^12 - V(type1, type2) / r^6);"
            "corr = (W(type1, type2) / rcut^12 - V(type1, type2) / rcut^6);"
            "ES = f/epsilon_r*q1*q2 * (1/r + krf * r^2 - crf);"
            "crf = 1 / rcut + krf * rcut^2;"
            "krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {self.epsilon_r};"
            "f = 138.935458;"
            f"rcut={self.nonbonded_cutoff.value_in_unit(nanometer)};"
        )
        # self._nb_force.setUsesPeriodicBoundaryConditions(True)
        # doesn't seem to exist for nb forces
        self._nb_force.addPerParticleParameter("type")
        self._nb_force.addPerParticleParameter("q")
        self._nb_force.setNonbondedMethod(
            mm.CustomNonbondedForce.CutoffPeriodic)
        self._nb_force.setCutoffDistance(self.nonbonded_cutoff
                                         .value_in_unit(nanometer))
        self._system.addForce(self._nb_force)
        self._forces_list.append(self._nb_force)

        for (_, type, charge, _) in filter(None, self._part_list):
            part_type_id = self.use_atom_type(type)
            self._nb_force.addParticle([part_type_id, charge])

        for (i, j) in filter(None, self._exclusion_list):
            self._nb_force.addExclusion(i, j)

        # add LJ parameters to the system
        C6 = []
        C12 = []
        # i,j => type index; t1,t2 => type names
        n = len(self._used_atom_types)
        for t1, i in self._used_atom_types.items():
            for t2, j in self._used_atom_types.items():
                nb_params = self._nb_types.get((t1, t2)) or\
                            self._nb_types.get((t2, t1)) or (0., 0.)
                V, W = nb_params
                c6 = 4 * W * (V ** 6)
                c12 = 4 * W * (V ** 12)
                C6.append(c6)
                C12.append(c12)
        self._nb_force.addTabulatedFunction(
            "V", mm.Discrete2DFunction(n, n, C6)
        )
        self._nb_force.addTabulatedFunction(
            "W", mm.Discrete2DFunction(n, n, C12)
        )
        self._nb_force_rebuild = False
        self._reinitialize = True

    def remove_force(self, force):
        for i, f in enumerate(self._forces_list):
            if f == force:
                del self._forces_list[i]
                self._system.removeForce(i)
                self._reinitialize = True
                return True
        return False

    def es_self_correction_add(self, i, j):
        _, _, q1, _ = self._part_list[i]
        _, _, q2, _ = self._part_list[j]
        qprod = q1 * q2
        if i == j:
            qprod *= 0.5
        if qprod != 0:
            self._es_self_correction_force.addBond(
                i, j, [qprod]
            )

    def build_es_self_correction_force(self):
        self._es_self_correction_force = mm.CustomBondForce(
            f"step(rcut-r) * ES;"
            f"ES = f*q_product/epsilon_r * (krf * r^2 - crf);"
            f"crf = 1 / rcut + krf * rcut^2;"
            f"krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {self.epsilon_r};"
            f"f = 138.935458;"
            f"rcut={self.nonbonded_cutoff.value_in_unit(nanometer)};"
        )
        # https://manual.gromacs.org/documentation/current/reference-manual/functions/nonbonded-interactions.html
        # see section Coulomb interaction with reaction field
        self._es_self_correction_force.addPerBondParameter("q_product")
        self._es_self_correction_force.setUsesPeriodicBoundaryConditions(True)
        self._system.addForce(self._es_self_correction_force)
        self._forces_list.append(self._es_self_correction_force)
        for i, (_, _, charge, _) in enumerate(filter(None, self._part_list)):
            if charge != 0:
                # self term in reaction field correction
                self.es_self_correction_add(i, i)
        for (i, j) in self._exclusion_list:
            self.es_self_correction_add(i, j)
        self._es_self_correction_force_rebuild = False
        self._reinitialize = True

    def add_exclusion(self, i, j):
        self._exclusion_list.append((i, j))
        if self.context_initialized:
            self._nb_force.addExclusion(i, j)
            self.es_self_correction_add(i, j)
            self._reinitialize = True
        return len(self._exclusion_list) - 1

    def get_exclusion_members(self, excl_id):
        return self._exclusion_list[excl_id]

    def remove_exclusion(self, excl_id):
        self._nb_force_rebuild = True
        self._exclusion_list[excl_id] = None

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
        self._nb_types[(type1, type2)] = (V, W)

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
        self.build_nb_force()
        self.build_es_self_correction_force()
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
        if self._es_self_correction_force_rebuild:
            self.remove_force(self._es_self_correction_force)
            self.build_es_self_correction_force()
        if self._nb_force_rebuild:
            self.remove_force(self._nb_force)
            self.build_nb_force()
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
