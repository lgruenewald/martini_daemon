#!/usr/bin/env python3

import openmm as mm
from openmm.unit import nanometer
import math
from collections import OrderedDict

class SysStar():

    _system: mm.System
    _context: mm.Context
    _nb_force: mm.Force
    _es_self_excl_force: mm.Force
    _es_except_force: mm.Force
    _lj_except_force: mm.Force
    _harmonic_angle_force: mm.Force
    _bond_force: mm.Force
    _proper_dihedral_force: mm.Force
    _improper_dihedral_force: mm.Force

    """All particles in the system
    indices are called part_id
    values are (part_name: string, part_type: string, charge: float,
    mass: float)"""
    _part_list: list[(str, str)]
    _used_atom_types: OrderedDict[str, int]
    context_initialized: bool

    def add_particle(self, part_name, part_type, charge, mass):
        """Adds a particle to the list, and returns its part_id"""
        defaults = self._atom_types.get(part_type)
        if defaults is None:
            raise ValueError("Attempt to add particle with unknown atom type"
                             f"Particle name: {part_name}. "
                             f"Particle type: {part_type}.")
        dcharge, dmass = defaults
        charge = charge or dcharge
        mass = mass or dmass
        self._part_list.append((part_name, part_type, charge, mass))
        self._system.addParticle(mass)
        part_type_id = self.use_atom_type(part_type)
        self._nb_force.addParticle([part_type_id, charge])
        return len(self._part_list) - 1

    def use_atom_type(self, part_type):
        """Add new atom types that might not exist in the first place.
        TODO: call this for every atom type that shows up in reaction products.

        It is not a problem to call this with atom types already present.

        Returns the atom type index
        """
        id = self._used_atom_types.get(part_type)
        if id is None:
            if self.context_initialized:
                raise Exception("use_atom_type call after context initialized")
            id = len(self._used_atom_types)
            self._used_atom_types[part_type] = id
        return id

    def len_particles(self):
        return len(self._part_list)

    """All bonds in the system
    indices are called bond_id
    values are (i: part_id, j: part_id, length: float, force: float)"""
    _harmonic_bond_list: list[(int, int, float, float)]

    def add_bond(self, part_id_i, part_id_j, length, force):
        """Adds a bond to the list, returns its bond_id"""
        if length is None or force is None:
            raise NotImplementedError("TODO length and strength for bonds must be specified for now")
        self._harmonic_bond_list.append((part_id_i, part_id_j, length, force))
        self._bond_force.addBond(part_id_i, part_id_j, length, force)
        return len(self._bond_list) - 1

    """All the angles in the system
    indices are called angle_id
    values are (i: part_id, j: part_id, k: part_id, theta: float, force: float)
    """
    _harmonic_angle_list: list[(int, int, int, float, float)]

    def add_angle(self, i, j, k, theta, force):
        """Adds an angle to the angle list, returns its angle_id"""
        self._harmonic_angle_list.append(
            (i, j, k, theta, force)
        )
        theta_rad = theta * math.pi / 180
        self._harmonic_angle_force.addAngle(i, j, k, theta_rad, force)
        return len(self._harmonic_angle_list) - 1

    """All the proper dihedrals in the system
    indices are called proper_dihedral_id
    values are (i: part_id, j: part_id, k: part_id, l: part_id,
     theta: float, force: float, multiplicity: int)"""
    _proper_dihedral_list: list[(int, int, int, int, float, float, int)]

    def add_proper_dihedral(self, i, j, k, l, theta, force, multiplicity):
        """Adds a dihedral to the list, returns its proper_dihedral_id"""
        self._proper_dihedral_list.append(
            (i, j, k, l, theta, force, multiplicity)
        )
        theta_rad = theta * math.pi / 180
        self._proper_dihedral_force.addTorsion(i, j, k, l, multiplicity,
                                               theta_rad, force)
        return len(self._proper_dihedral_list) - 1

    """All the exclusions in the system
    indices are caled excl_id
    values are (i: part_id, j: part_id)
    """
    _exclusion_list: list[(int, int)]

    def add_exclusion(self, i, j):
        self._exclusion_list.append((i, j))
        return len(self._exclusion_list) - 1

    """All the constraints in the system
    indices are called constraint_id
    values are (i: part_id, j: part_id, length: float)
    """
    _constraint_list: list[(int, int, float)]

    def add_constraint(self, i, j, length):
        """Adds a constraint to the list, returns its constraint_id"""
        self._constraint_list.append((i, j, length))
        return len(self._constraint_list) - 1

    """All the improper dihedrals in the system
    indices are called improper_dihedral_id
    values are (i: part_id, j: part_id, k: part_id, l: part_id,
     theta: float, force: float)
    """
    _improper_dihedral_list: list[(int, int, int, int, float, float)]

    def add_improper_dihedral(self, i, j, k, l, theta, force):
        """Adds a dihedral to the angle list, returns its improper_dihedral_id
        """
        self._improper_dihedral_list.append(
            (i, j, k, l, theta, force)
        )
        theta = theta - 360 if theta > 180 else theta
        theta_rad = theta * math.pi / 180
        self._improper_dihedral_force.addTorsion(i, j, k, l,
                                                 (theta_rad, force))
        return len(self._improper_dihedral_list) - 1

    """Atom types to look up default charges and default masses
    values are atom_type: string, (mass: float, charge: float)
    """
    _atom_types: dict[str, (float, float)]

    def add_atom_type(self, type, mass, charge):
        self._atom_types[type] = (mass, charge)

#    _bond_types: list  # TODO
#    _angle_types: list  # TODO
#    _dihedral_types: list  # TODO

    """Non bonded parameters are for building the C6/C12 table
    values are (type1: string, type2: string), (V: float, W: float)
    """
    _nb_types: dict[(str, str), (float, float)]

    def add_nb_type(self, type1, type2, V, W):
        self._nb_types[(type1, type2)] = (V, W)

    def __init__(self):
        self._system = mm.System()
        self._part_list = []
        self._used_atom_types = {}
        self._harmonic_bond_list = []
        self._harmonic_angle_list = []
        self._proper_dihedral_list = []
        self._exclusion_list = []
        self._constraint_list = []
        self._improper_dihedral_list = []
        self._atom_types = OrderedDict()
        self.context_initialized = False
        self._nb_types = {}
        self.epsilon_r = 15.0  # TODO unhardcode
        self.nonbonded_cutoff = 1.1 * nanometer
        self.nonbonded_cutoff
        # misc forces
        self._system.addForce(mm.CMMotionRemover())
        # Non bonded force
        self._nb_force = mm.CustomNonbondedForce(
            "step(rcut-r)*(LJ - corr + ES);"
            "LJ = (C12(type1, type2) / r^12 - C6(type1, type2) / r^6);"
            "corr = (C12(type1, type2) / rcut^12 - C6(type1, type2) / rcut^6);"
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

        # TODO tabulated nb force params

        self._bond_force = mm.HarmonicBondForce()
        self._system.addForce(self._bond_force)

        self._harmonic_angle_force = mm.HarmonicAngleForce()
        self._system.addForce(self._harmonic_angle_force)

        self._proper_dihedral_force = mm.PeriodicTorsionForce()
        self._system.addForce(self._proper_dihedral_force)

        self._improper_dihedral_force = mm.CustomTorsionForce(
            "0.5*k*(thetap-theta0)^2; thetap = "
            "step(-(theta-theta0+pi))*2*pi+theta+step(theta-theta0-pi)*(-2*pi)"
            f"; pi = {math.pi:.14f}"
        )
        self._improper_dihedral_force.addPerTorsionParameter("theta0")
        self._improper_dihedral_force.addPerTorsionParameter("k")
        self._system.addForce(self._improper_dihedral_force)

        # TODO CMAPTorsionForce

    def build_context(self, integrator, periodicBoxVectors):
        """After build_context the following things should not happen:
        - Adding new *_types to _atom_types, _nb_types etc.
        - Using new atom types (through new particles with new
        types or calling use_atom_type) -- call use_atom_type for every type
        that will eventually be needed first
        """
        # Finish setup
        self.context_initialized = True
        # add LJ parameters to the system
        Vs = []
        Ws = []
        # i,j => type index; t1,t2 => type names
        n = len(self._used_atom_types)
        for i, t1 in self._used_atom_types.items():
            for j, t2 in self._used_atom_types.items():
                nb_params = self._nb_types[(t1, t2)] or\
                            self._nb_types[(t2, t1)]
                if nb_params is None:
                    raise ValueError(f"Couldn't find LJ params for {t1}; {t2}")
                V, W = nb_params
                Vs.append(V)
                Ws.append(W)
        self._nb_force.addTabulatedFunction(
            "C6", mm.Discrete2DFunction(n, n, Vs)
        )
        self._nb_force.addTabulatedFunction(
            "C12", mm.Discrete2DFunction(n, n, Ws)
        )

        # Build context
        self._context = mm.Context(self.system, integrator)
        self._context.setPeriodicBoxVectors(periodicBoxVectors)
        # === API TODOs to parallel openmm.app's Simulation ===
        # TODO positions, initial velocity, energy minimizations, couplings
        # TODO stepping the simulation forward
        # TODO reporters, getState checkpoints, ...

    def dump(self):
        print("Particles:", self._part_list)
        print("Bonds:", self._harmonic_bond_list)
        print("Angles:", self._harmonic_angle_list)
        print("Dihedrals:", self._proper_dihedral_list)
        print("Impropers:", self._improper_dihedral_list)
        print("Exclusions:", self._exclusion_list)
        print("Constraints:", self._constraint_list)
