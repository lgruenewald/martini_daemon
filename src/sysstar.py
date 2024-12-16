#!/usr/bin/env python3

from __future__ import annotations
import openmm as mm
import openmm.app as mmapp
from openmm.unit import nanometer, picosecond, md_unit_system
import math
import utils
from collections import OrderedDict
import numpy as np
from dataclasses import dataclass


class Force():
    _force_obj: mm.Force
    _sysstar: SysStar
    _list: list
    _rebuild: bool
    # if rebuild is set to True, it means that the force object _force_obj is
    # no longer valid, or does not exist. If rebuild is false, as much effort
    # should be done to keep _force_obj updated as possible
    # rebuild is True for example before force_obj is built in the first place,
    # or when removing elements from the bond

    def __init__(self, sysstar):
        self._list = []
        self._sysstar = sysstar
        self._rebuild = True
        self._force_obj = None

    def _build(self):
        pass

    def add(self, *params):
        pass

    def get_members(self, i):
        pass

    def update_params(self, i, *params):
        pass

    def remove(self, i):
        self._list[i] = None
        self._rebuild = True

    def build(self):
        if self._rebuild:
            self.destroy()
            self._build()
            self._rebuild = False
            self._sysstar._forces_list.append(self._force_obj)
            self._sysstar._system.addForce(self._force_obj)
            self._sysstar._reinitialize = True

    def destroy(self):
        if self._force_obj is None:
            return False
        for i, f in enumerate(self._sysstar._forces_list):
            if f == self._force_obj:
                del self._sysstar._forces_list[i]
                self._sysstar._system.removeForce(i)
                self._sysstar._reinitialize = True
                return True
        return False

    def _interaction(self):
        return Interaction(self, len(self._list) - 1)


@dataclass
class Interaction():
    _force: Force
    _index: int

    def get_members(self):
        return self._force.get_members(self._index)

    def update_params(self, *params):
        return self._force.update_params(self._index, params)

    def remove(self):
        self._force.remove(self._index)


class HarmonicBond(Force):
    """All bonds in the system
    indices are called bond_id
    values are (i: part_id, j: part_id, length: float, force: float)"""

    def _build(self):
        self._force_obj = mm.HarmonicBondForce()
        for (i, j, length, force) in filter(None, self._list):
            self._force_obj.addBond(i, j, length, force)

    def add(self, i, j, length, force):
        if length is None or force is None:
            raise NotImplementedError("Bond length and force are mandatory.")
        self._list.append((i, j, length, force))
        if not self._rebuild:
            self._force_obj.addBond(i, j, length, force)
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, _, _ = self._list[id]
        return [i, j]

    def update_params(self, id, length, force):
        raise NotImplementedError


class HarmonicAngle(Force):
    """All the angles in the system
    indices are called angle_id
    values are (i: part_id, j: part_id, k: part_id, theta: float, force: float)
    """

    def _build(self):
        self._force_obj = mm.HarmonicAngleForce()
        for (i, j, k, theta, force) in filter(None, self._list):
            self._force_obj.addAngle(i, j, k, theta, force)

    def add(self, i, j, k, theta, force):
        self._list.append((i, j, k, theta, force))
        if not self._rebuild:
            self._force_obj.addAngle(i, j, k, theta, force)
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, _, _ = self._list[id]
        return [i, j, k]

    def update_params(self, id, theta, force):
        raise NotImplementedError


class G96Angle(Force):
    def _build(self):
        self._force_obj = mm.CustomAngleForce(
            "0.5 * k * (cos(theta) - cos(theta0))^2"
        )
        self._force_obj.addPerAngleParameter("theta0")
        self._force_obj.addPerAngleParameter("k")
        for (i, j, k, theta, force) in filter(None, self._list):
            self._force_obj.addAngle(i, j, k, (theta, force))

    def add(self, i, j, k, theta, force):
        self._list.append((i, j, k, theta, force))
        if not self._rebuild:
            self._force_obj.addAngle(i, j, k, (theta, force))
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, _, _ = self._list[id]
        return [i, j, k]

    def update_params(self, id, theta, force):
        raise NotImplementedError


class RestrictedAngle(Force):
    def _build(self):
        self._force_obj = mm.CustomAngleForce(
            "0.5*k*(cos(theta)-cos(theta0))^2/sin(theta)^2"
        )
        self._force_obj.addPerAngleParameter("theta0")
        self._force_obj.addPerAngleParameter("k")
        for (i, j, k, theta, force) in filter(None, self._list):
            self._force_obj.addAngle(i, j, k, (theta, force))

    def add(self, i, j, k, theta, force):
        self._list.append((i, j, k, theta, force))
        if not self._rebuild:
            self._force_obj.addAngle(i, j, k, (theta, force))
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, _, _ = self._list[id]
        return [i, j, k]

    def update_params(self, id, theta, force):
        raise NotImplementedError


class ProperDihedral(Force):
    """All the proper dihedrals in the system
    indices are called proper_dihedral_id
    values are (i: part_id, j: part_id, k: part_id, l: part_id,
     theta: float, force: float, multiplicity: int)"""

    def _build(self):
        self._force_obj = mm.PeriodicTorsionForce()
        for (i, j, k, l, theta, force, mult) in filter(None, self._list):
            self._force_obj.addTorsion(i, j, k, l, mult, theta, force)

    def add(self, i, j, k, l, theta, force, mult):
        self._list.append((i, j, k, l, theta, force, mult))
        if not self._rebuild:
            self._force_obj.addTorsion(i, j, k, l, mult, theta, force)
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, l, _, _, _ = self._list[id]
        return [i, j, k, l]

    def update_params(self, id, theta, force, mult):
        raise NotImplementedError


class ImproperDihedral(Force):
    """All the improper dihedrals in the system
    indices are called improper_dihedral_id
    values are (i: part_id, j: part_id, k: part_id, l: part_id,
     theta: float, force: float)
    """

    def _build(self):
        self._force_obj = mm.CustomTorsionForce(
            "0.5*k*(thetap-theta0)^2;"
            "thetap = step(-plus)*2*pi+theta+step(minus)*(-2*pi);"
            "plus=theta+pi-theta0;"
            "minus=theta-pi-theta0;"
            f"pi = {math.pi:.14f}"
        )
        self._force_obj.addPerTorsionParameter("theta0")
        self._force_obj.addPerTorsionParameter("k")
        for (i, j, k, l, theta, force) in filter(None, self._list):
            self._force_obj.addTorsion(i, j, k, l, (theta, force))

    def add(self, i, j, k, l, theta, force):
        self._list.append((i, j, k, l, theta, force))
        if not self._rebuild:
            self._force_obj.addTorsion(i, j, k, l, (theta, force))
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, l, _, _ = self._list[id]
        return [i, j, k, l]

    def update_params(self, id, theta, force):
        raise NotImplementedError


class CombinedBendingTorsion(Force):
    def _build(self):
        self._force_obj = mm.CustomCompoundBondForce(
            4,
            "k*sintheta0^3*sintheta1^3*(a0 + a1*cosphi + a2*cosphi^2 + a3*cosphi^3 + a4*cosphi^4); "
            "sintheta0 = sin(angle(p1, p2, p3));"
            "sintheta1 = sin(angle(p2, p3, p4));"
            "cosphi = cos(dihedral(p1, p2, p3, p4));",
        )
        self._force_obj.addPerBondParameter("k")
        self._force_obj.addPerBondParameter("a0")
        self._force_obj.addPerBondParameter("a1")
        self._force_obj.addPerBondParameter("a2")
        self._force_obj.addPerBondParameter("a3")
        self._force_obj.addPerBondParameter("a4")
        for (i, j, k, l, force, a0, a1, a2, a3, a4) in filter(None, self._list):
            self._force_obj.addBond((i, j, k, l), (force, a0, a1, a2, a3, a4))

    def add(self, i, j, k, l, force, a0, a1, a2, a3, a4):
        self._list.append((i, j, k, l, force, a0, a1, a2, a3, a4))
        if not self._rebuild:
            self._force_obj.addBond((i, j, k, l), (force, a0, a1, a2, a3, a4))
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, l, _, _, _, _, _, _ = self._list[id]
        return [i, j, k, l]

    def update_params(self, id, theta, force):
        raise NotImplementedError


class RBTorsion(Force):
    def _build(self):
        self._force_obj = mm.RBTorsionForce()
        for (i, j, k, l, c0, c1, c2, c3, c4, c5) in filter(None, self._list):
            self._force_obj.addTorsion(i, j, k, l, c0, c1, c2, c3, c4, c5)

    def add(self, i, j, k, l, c0, c1, c2, c3, c4, c5):
        self._list.append((i, j, k, l, c0, c1, c2, c3, c4, c5))
        if not self._rebuild:
            self._force_obj.addTorsion(i, j, k, l, c0, c1, c2, c3, c4, c5)
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, l, _, _, _ = self._list[id]
        return [i, j, k, l]

    def update_params(self, id, theta, force, mult):
        raise NotImplementedError


class VSiteWeighedAverage(Force):
    # Force because it implements the same interface, but it reimplements
    # everything, the default helpers don't suit it well (no _force_obj).

    def _make_vsite(self, vid, members, weights):
        n = len(members)
        match n:
            case 1:
                vsite = mm.TwoParticleAverageSite(
                    members[0], members[0], 1.0, 0.0
                )
            case 2:
                vsite = mm.TwoParticleAverageSite(
                    members[0], members[1], weights[0], weights[1]
                )
            case 3:
                vsite = mm.ThreeParticleAverageSite(
                    members[0], members[1], members[2],
                    weights[0], weights[1], weights[2]
                )
            case _:
                vsite = mm.LocalCoordinatesSite(
                    members,  # particles
                    weights,  # origin weights
                    [0.0] * n,  # x direction weight
                    [0.0] * n,  # y direction weight
                    [0.0, 0.0, 0.0]  # coordinates
                )
        self._sysstar._system.setVirtualSite(vid, vsite)

    def _build(self):
        for (vid, members, weights) in self._list:
            self._make_vsite(vid, members, weights)

    def add(self, members, weights):
        vid, *members = members
        self._list.append((vid, members, weights))
        if not self._rebuild:
            self._make_vsite(vid, members, weights)
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, i):
        vid, members, _ = self._list[i]
        return [vid, *members]

    def update_params(self, i, *params):
        raise NotImplementedError

    def remove(self, i):
        raise NotImplementedError

    def build(self):
        # only should get called once when building it initially
        if self._rebuild:
            self._build()
            self._rebuild = False
            self._sysstar._reinitialize = True

    def destroy(self):
        # no _force_obj, so it shold never be called like this
        raise NotImplementedError

    def _interaction(self):
        return Interaction(self, len(self._list) - 1)


class VSite3fad(Force):
    # Force because it implements the same interface, but it reimplements
    # everything, the default helpers don't suit it well (no _force_obj).

    def _make_vsite(self, vid, i, j, k, theta, d):
        dcos = d * math.cos(theta)
        dsin = d * math.sin(theta)
        vsite = mm.LocalCoordinatesSite(
            [i, j, k],  # particles
            [1.0, 0.0, 0.0],  # origin weights
            [-0.5, 0.5, 0.0],  # x direction weight
            [0.0, -0.5, 0.5],  # y direction weight
            [dcos, dsin, 0.0]  # coordinates
        )
        self._sysstar._system.setVirtualSite(vid, vsite)

    def _build(self):
        for (vid, i, j, k, theta, d) in self._list:
            self._make_vsite(vid, i, j, k, theta, d)

    def add(self, members, params):
        vid, i, j, k = members
        theta, d = params
        self._list.append((vid, i, j, k, theta, d))
        if not self._rebuild:
            self._make_vsite(vid, i, j, k, theta, d)
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, i):
        vid, i, j, k, _, _ = self._list[i]
        return [vid, i, j, k]

    def update_params(self, i, *params):
        raise NotImplementedError

    def remove(self, i):
        raise NotImplementedError

    def build(self):
        # only should get called once when building it initially
        if self._rebuild:
            self._build()
            self._rebuild = False
            self._sysstar._reinitialize = True

    def destroy(self):
        # no _force_obj, so it shold never be called like this
        raise NotImplementedError

    def _interaction(self):
        return Interaction(self, len(self._list) - 1)


class VSite3out(Force):
    # Force because it implements the same interface, but it reimplements
    # everything, the default helpers don't suit it well (no _force_obj).

    def _make_vsite(self, vid, i, j, k, a, b, c):
        vsite = mm.OutOfPlaneSite(
            i, j, k, a, b, c
        )
        self._sysstar._system.setVirtualSite(vid, vsite)

    def _build(self):
        for (vid, i, j, k, a, b, c) in self._list:
            self._make_vsite(vid, i, j, k, a, b, c)

    def add(self, members, params):
        vid, i, j, k = members
        a, b, c = params
        self._list.append((vid, i, j, k, a, b, c))
        if not self._rebuild:
            self._make_vsite(vid, i, j, k, a, b, c)
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, i):
        vid, i, j, k, _, _ = self._list[i]
        return [vid, i, j, k]

    def update_params(self, i, *params):
        raise NotImplementedError

    def remove(self, i):
        raise NotImplementedError

    def build(self):
        # only should get called once when building it initially
        if self._rebuild:
            self._build()
            self._rebuild = False
            self._sysstar._reinitialize = True

    def destroy(self):
        # no _force_obj, so it shold never be called like this
        raise NotImplementedError

    def _interaction(self):
        return Interaction(self, len(self._list) - 1)


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

    def __init__(self):
        self._part_list = []
        self._used_atom_types = {}
        self._exclusion_list = []
        self._constraint_list = []
        self._atom_types = OrderedDict()
        self.context_initialized = False
        self._reinitialize = True
        self._nb_types = {}
        self._forces_list = []
        self.epsilon_r = 15.0  # TODO unhardcode
        self.nonbonded_cutoff = 1.1 * nanometer
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
            self._system.addConstraint(i, j, length)
            self._reinitialize = True
        return len(self._constraint_list) - 1

    def get_constraint_members(self, constraint_id):
        i, j, _ = self._constraint_list[constraint_id]
        return (i, j)

    def remove_constraint(self, constraint_id):
        raise NotImplementedError("Constraints cannot be safely removed yet")

    def add_atom_type(self, type, charge, mass):
        if self.context_initialized:
            raise ValueError("Cannot do this after context is initialized")
        self._atom_types[type] = (charge, mass)

    def add_nb_type(self, type1, type2, V, W):
        if self.context_initialized:
            raise ValueError("Cannot do this after context is initialized")
        self._nb_types[(type1, type2)] = (V, W)

    def build_context(self, integrator, periodicBoxVectors):
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
        self._context = mm.Context(self._system, integrator)
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
        box = state.getPeriodicBoxVectors()[0].x
        pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)
        pos = np.remainder(pos, box)
        return pos, box

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
        self._context.applyConstraints(tol=1e-5)

    _xtc: mmapp.XTCFile = None

    def set_xtc_path(self, path):
        if not self.context_initialized:
            raise Exception("Initialize the context first")
        utils.backup_try(path)
        mmtopol = mmapp.Topology()
        mmtopol._numAtoms = self.len_particles()
        mmtopol._periodicBoxVectors = self._periodic_box
        timestep = self._integrator.getStepSize()
        self._xtc = mmapp.XTCFile(path, mmtopol, timestep)

    def write_gro(self, path):
        if not self.context_initialized:
            raise Exception("Initialize the context first")
        utils.backup_try(path)
        state = self._context.getState(positions=True, velocities=True)
        box = state.getPeriodicBoxVectors()[0].x
        pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)
        pos = np.remainder(pos, box)
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
        utils.backup_try("sys.dump")
        with open("sys.dump", "w") as file:
            print("==== SysStar Dump ====", file=file)
            print("Particles:", self._part_list, file=file)
            for force in self.modular_forces:
                print(force.__class__.__name__, force._list, file=file)
            print("Exclusions:", self._exclusion_list, file=file)
            print("Constraints:", self._constraint_list, file=file)
