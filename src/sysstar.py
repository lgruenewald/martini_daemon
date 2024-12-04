#!/usr/bin/env python3

import openmm as mm
import openmm.app as mmapp
from openmm.unit import nanometer, picosecond
import math
import utils
from collections import OrderedDict

class SysStar():

    _system: mm.System
    _context: mm.Context
    _nb_force: mm.Force = None
    _nb_force_rebuild = False
    _es_self_correction_force: mm.Force = None
    _es_self_correction_force_rebuild = False
    _harmonic_angle_force: mm.Force = None
    _harmonic_angle_rebuild = False
    _bond_force: mm.Force = None
    _bond_rebuild = False
    _proper_dihedral_force: mm.Force = None
    _proper_dihedral_rebuild = False
    _improper_dihedral_force: mm.Force = None
    _improper_dihedral_rebuild = False
    _forces_list: list[mm.Force]  # to keep track of indices

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
        i = len(self._part_list) - 1
        if charge != 0:
            # self term in reaction field correction
            # TODO map this part_id -> index so charges can be changed reliably
            self._es_self_correction_force.addBond(i, i, [0.5 * charge ** 2])
        return i

    def get_particle_name_type(self, i):
        name, type, _, _ = self._part_list[i]
        return name, type

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
        return len(self._harmonic_bond_list) - 1

    def get_bond_members(self, bond_id):
        i, j, _, _ = self._harmonic_bond_list[bond_id]
        return (i, j)

    def remove_bond(self, bond_id):
        self._bond_rebuild = True
        self._harmonic_bond_list[bond_id] = None

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

    def get_angle_members(self, angle_id):
        i, j, k, _, _ = self._harmonic_angle_list[angle_id]
        return (i, j, k)

    def remove_angle(self, angle_id):
        self._harmonic_angle_rebuild = True
        self._harmonic_angle_list[angle_id] = None

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

    def get_proper_dihedral_members(self, dih_id):
        i, j, k, l, _, _, _ = self._proper_dihedral_list[dih_id]
        return (i, j, k, l)

    def remove_proper_dihedral(self, dih_id):
        self._proper_dihedral_rebuild = True
        self._proper_dihedral_list[dih_id] = None

    """All the exclusions in the system
    indices are caled excl_id
    values are (i: part_id, j: part_id)
    """
    _exclusion_list: list[(int, int)]

    def add_exclusion(self, i, j):
        self._exclusion_list.append((i, j))
        self._nb_force.addExclusion(i, j)
        return len(self._exclusion_list) - 1

    def get_exclusion_members(self, excl_id):
        return self._exclusion_list[excl_id]

    def remove_exclusion(self, excl_id):
        self._nb_force_rebuild = True
        self._exclusion_list[excl_id] = None

    """All the constraints in the system
    indices are called constraint_id
    values are (i: part_id, j: part_id, length: float)
    """
    _constraint_list: list[(int, int, float)]

    def add_constraint(self, i, j, length):
        """Adds a constraint to the list, returns its constraint_id"""
        self._constraint_list.append((i, j, length))
        self._system.addConstraint(i, j, length)
        return len(self._constraint_list) - 1

    def get_constraint_members(self, constraint_id):
        i, j, _ = self._constraint_list[constraint_id]
        return (i, j)

    def remove_constraint(self, constraint_id):
        self._system.removeConstraint(constraint_id)
        # TODO does this mangle the IDs of constraints that come after?
        self._constraint_list[constraint_id] = None

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

    def get_improper_members(self, dih_id):
        i, j, k, l, _, _ = self._improper_dihedral_list[dih_id]
        return (i, j, k, l)

    def remove_improper(self, dih_id):
        self._improper_dihedral_rebuild = True
        self._improper_dihedral_list[dih_id] = None

    """Atom types to look up default charges and default masses
    values are atom_type: string, (mass: float, charge: float)
    """
    _atom_types: dict[str, (float, float)]

    def add_atom_type(self, type, charge, mass):
        self._atom_types[type] = (charge, mass)

#    _bond_types: list  # TODO
#    _angle_types: list  # TODO
#    _dihedral_types: list  # TODO

    """Non bonded parameters are for building the C6/C12 table
    values are (type1: string, type2: string), (V: float, W: float)
    """
    _nb_types: dict[(str, str), (float, float)]

    def add_nb_type(self, type1, type2, V, W):
        self._nb_types[(type1, type2)] = (V, W)

    _context: mm.Context
    _integrator: mm.Integrator

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

        # TODO tabulated nb force params

        self._bond_force = mm.HarmonicBondForce()
        self._system.addForce(self._bond_force)
        self._forces_list.append(self._bond_force)

        self._harmonic_angle_force = mm.HarmonicAngleForce()
        self._system.addForce(self._harmonic_angle_force)
        self._forces_list.append(self._harmonic_angle_force)

        self._proper_dihedral_force = mm.PeriodicTorsionForce()
        self._system.addForce(self._proper_dihedral_force)
        self._forces_list.append(self._proper_dihedral_force)

        self._improper_dihedral_force = mm.CustomTorsionForce(
            "0.5*k*(thetap-theta0)^2; thetap = "
            "step(-(theta-theta0+pi))*2*pi+theta+step(theta-theta0-pi)*(-2*pi)"
            f"; pi = {math.pi:.14f}"
        )
        self._improper_dihedral_force.addPerTorsionParameter("theta0")
        self._improper_dihedral_force.addPerTorsionParameter("k")
        self._system.addForce(self._improper_dihedral_force)
        self._forces.append(self._improper_dihedral_force)

        self._es_self_correction_force = mm.CustomBondForce(
            f"step(rcut-r) * ES;"
            f"ES = f*q_product/epsilon_r * (krf * r^2 - crf);"
            f"crf = 1 / rcut + krf * rcut^2;"
            f"krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {self.epsilon_r};"
            f"f = 138.935458;"
            f"rcut={self.nonbonded_cutoff.value_in_unit(nanometer)};"
        )
        self._es_self_correction_force.addPerBondParameter("q_product")
        self._system.addForce(self._es_self_correction_force)
        self._forces.append(self._es_self_correction_force)

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
        self._integrator = integrator
        # add LJ parameters to the system
        C6 = []
        C12 = []
        # TODO fixme why are they sometimes called V,W or C6,C12
        # i,j => type index; t1,t2 => type names
        n = len(self._used_atom_types)
        for t1, i in self._used_atom_types.items():
            for t2, j in self._used_atom_types.items():
                nb_params = self._nb_types.get((t1, t2)) or\
                            self._nb_types.get((t2, t1))
                if nb_params is None:
                    raise ValueError(f"Couldn't find LJ params for {t1}; {t2}")
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

        # Build context
        self._periodic_box = periodicBoxVectors
        self._system.setDefaultPeriodicBoxVectors(*periodicBoxVectors)
        self._context = mm.Context(self._system, integrator)
        self._context.setPeriodicBoxVectors(*periodicBoxVectors)

    def reinitialize(self):
        if not self.context_initialized:
            raise Exception("Initialize the context first")
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
        if self.context_initialized:
            raise Exception("Add forces before initializing the context")
        self._system.addForce(force)

    def do_steps(self, steps):
        if not self.context_initialized:
            raise Exception("Initialize the context first")
        self._integrator.step(steps)
        if self._xtc is not None:
            self._xtc.interval = steps
            pos = self.get_state().getPositions()
            self._xtc.writeModel(pos)

    def get_state(self):
        if not self.context_initialized:
            raise Exception("Initialize the context first")
        return self._context.getState(positions=True, velocities=True,
                                      forces=True, energy=True,
                                      enforcePeriodicBox=True)

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
        state = self.get_state()
        pos = state.getPositions()
        vel = state.getVelocities()
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
                           f"{atom_index:5}{cpos.x:8.3f}{cpos.y:8.3f}"
                           f"{cpos.z:8.3f}{cvel.x:8.4f}{cvel.y:8.4f}"
                           f"{cvel.z:8.4f}\n")
            v1, v2, v3 = (v.value_in_unit(nanometer)
                          for v in state.getPeriodicBoxVectors())
            file.write(
                f"{v1[0]:.4f} {v2[1]:.4f} {v3[2]:.4f} {v1[1]:.4f} {v1[2]:.4f} "
                f"{v2[0]:.4f} {v2[2]:.4f} {v3[0]:.4f} {v3[1]:.4f}\n"
            )

    def dump(self):
        print("==== SysStar Dump ====")
        print("Particles:", self._part_list)
        print("Bonds:", self._harmonic_bond_list)
        print("Angles:", self._harmonic_angle_list)
        print("Dihedrals:", self._proper_dihedral_list)
        print("Impropers:", self._improper_dihedral_list)
        print("Exclusions:", self._exclusion_list)
        print("Constraints:", self._constraint_list)
        print("Used atom types: ", self._used_atom_types)
        Vs = []
        Ws = []
        # i,j => type index; t1,t2 => type names
        for t1, i in self._used_atom_types.items():
            for t2, j in self._used_atom_types.items():
                nb_params = self._nb_types.get((t1, t2)) or\
                            self._nb_types.get((t2, t1))
                if nb_params is None:
                    raise ValueError(f"Couldn't find LJ params for {t1}; {t2}")
                V, W = nb_params
                Vs.append(V)
                Ws.append(W)
        print("nb params V", Vs)
        print("nb params W", Ws)
