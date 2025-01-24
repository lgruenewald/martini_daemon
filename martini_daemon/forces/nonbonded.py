from __future__ import annotations
import openmm as mm
from openmm.unit import nanometer
from .force import Force, Interaction
from collections import OrderedDict


class NonBonded(Force):
    """
    S* modular_forces implementation of non bonded forces.
    Different to other Force objects because of the interplay between the
    NB force, self exclusion force, exclusions and atomtypes.
    """

    """Non bonded parameters are for building the C6/C12 table
    values are (type1: string, type2: string), (V: float, W: float)
    """
    _nb_types: dict[(str, str), (float, float)]
    _used_atom_types: OrderedDict[str, int]
    _exclusions: ExclusionHelper
    _rebuild: bool
    _force_obj: mm.Force

    def __init__(self, sysstar):
        self._sysstar = sysstar
        self._rebuild = True
        self._force_obj = None
        self._used_atom_types = {}
        self._nb_types = {}
        self._exclusions = ExclusionHelper(self, sysstar)

    def add(self, *params):
        """Particle types are listed in S*, not here"""
        raise NotImplementedError

    def get_members(self, i):
        """Non bonded force, so invalid"""
        raise NotImplementedError

    def remove(self, i):
        """Non bonded force, so invalid"""
        raise NotImplementedError

    def update_params(self, i, part_type, charge, charge_changed):
        part_type_id = self.use_atom_type(part_type)
        self._force_obj.setPerParticleParameters(i, [part_type_id, charge])
        self._sysstar.pairs._rebuild = True
        if charge_changed:
            self._es_force._rebuild = True

    def _build(self):
        self._force_obj = mm.CustomNonbondedForce(
            "step(rcut-r)*(LJ - corr + ES);"
            "LJ = (W(type1, type2) / r^12 - V(type1, type2) / r^6);"
            "corr = (W(type1, type2) / rcut^12 - V(type1, type2) / rcut^6);"
            "ES = f/epsilon_r*q1*q2 * (1/r + krf * r^2 - crf);"
            "crf = 1 / rcut + krf * rcut^2;"
            "krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {self._sysstar.epsilon_r};"
            "f = 138.935458;"
            f"rcut={self._sysstar.nonbonded_cutoff.value_in_unit(nanometer)};"
        )
        self._force_obj.addPerParticleParameter("type")
        self._force_obj.addPerParticleParameter("q")
        self._force_obj.setNonbondedMethod(
            mm.CustomNonbondedForce.CutoffPeriodic
        )
        self._force_obj.setCutoffDistance(
            self._sysstar.nonbonded_cutoff.value_in_unit(nanometer)
        )

        for (_, type, charge, _) in filter(None, self._sysstar._part_list):
            part_type_id = self.use_atom_type(type)
            self._force_obj.addParticle([part_type_id, charge])

        for (i, j) in filter(None, self._exclusions._list):
            self._force_obj.addExclusion(i, j)

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
        self._force_obj.addTabulatedFunction(
            "V", mm.Discrete2DFunction(n, n, C6)
        )
        self._force_obj.addTabulatedFunction(
            "W", mm.Discrete2DFunction(n, n, C12)
        )

    def build(self):
        if self._rebuild:
            self.destroy()
            self._build()
            # self._force_obj.setUsesPeriodicBoundaryConditions(True)
            # not applicable for non bonded, they use setNonBondedMethod
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
        raise NotImplementedError

    def get_exclusion_helper(self):
        return self._exclusions

    def use_atom_type(self, part_type):
        """Add new atom types for LJ.

        It is not a problem to call this with atom types already present.

        Returns the atom type index
        """
        id = self._used_atom_types.get(part_type)
        if id is None:
            id = len(self._used_atom_types)
            self._used_atom_types[part_type] = id
            self._rebuild = True
        return id


class ExclusionHelper(Force):
    """
    Just a wrapper around handling exclusions using Interaction
    handled by NonBonded

    Also handles the electrostatic self correction force

    Different to other Force objects because it is partly managed by NonBonded.

    See:
    https://manual.gromacs.org/documentation/current/reference-manual/functions/nonbonded-interactions.html
    """
    # TODO is a dict based approach better, then duplicate additions of
    # exclusions could return the same one
    # does returning an existing exclusion create any problems? what if that
    # exclusion gets removed in a reaction?
    _list: list[(int, int)]
    _nb: NonBonded
    _rebuild: bool
    _force_obj: mm.Force

    def __init__(self, nb, sysstar):
        self._list = []
        self._nb = nb
        self._sysstar = sysstar
        self._rebuild = True
        self._force_obj = None

    def es_self_correction_add(self, i, j):
        _, _, q1, _ = self._sysstar._part_list[i]
        _, _, q2, _ = self._sysstar._part_list[j]
        qprod = q1 * q2
        if i == j:
            qprod *= 0.5
        if qprod != 0:
            self._force_obj.addBond(
                i, j, [qprod]
            )

    def _build(self):
        self._force_obj = mm.CustomBondForce(
            f"step(rcut-r) * ES;"
            f"ES = f*q_product/epsilon_r * (krf * r^2 - crf);"
            f"crf = 1 / rcut + krf * rcut^2;"
            f"krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {self._sysstar.epsilon_r};"
            f"f = 138.935458;"
            f"rcut={self._sysstar.nonbonded_cutoff.value_in_unit(nanometer)};"
        )
        self._force_obj.addPerBondParameter("q_product")
        self._force_obj.setUsesPeriodicBoundaryConditions(True)
        for i, (_, _, charge, _) in \
                enumerate(filter(None, self._sysstar._part_list)):
            if charge != 0:
                # self term in reaction field correction
                self.es_self_correction_add(i, i)
        for (i, j) in self._list:
            self.es_self_correction_add(i, j)

    def add(self, i, j):
        self._list.append((i, j))
        if not self._rebuild:
            self.es_self_correction_add(i, j)
            self._sysstar._reinitialize = True
        if not self._nb._rebuild:
            self._nb._force_obj.addExclusion(i, j)
            self._sysstar._reinitialize = True

        return self._interaction()

    def get_members(self, i):
        return self._list[i]

    def update_params(self, i, *params):
        raise NotImplementedError

    def remove(self, i):
        self._list[i] = None
        self._rebuild = True
        self._nb._rebuild = True

    def build(self):
        # dependency on NB being built first
        if self._nb._rebuild:
            self._nb.build()
        if self._rebuild:
            self.destroy()
            self._build()
            self._force_obj.setUsesPeriodicBoundaryConditions(True)
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
