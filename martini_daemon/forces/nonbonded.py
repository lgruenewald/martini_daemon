from __future__ import annotations
import openmm as mm
from .force import Force, Interaction
from collections import OrderedDict


class NonBonded(Force):
    """
    S* modular_forces implementation of non bonded forces.
    Different to other Force objects because of the interplay between the
    NB force, self exclusion force, exclusions and atomtypes.
    """

    def __init__(self, epsilon_r=15.0, cutoff_nm=1.1):
        self.epsilon_r = epsilon_r
        self.cutoff_nm = cutoff_nm
        self._rebuild: bool = True
        self._force_obj: None | mm.Force = None
        self._used_atom_types: OrderedDict[str, int] = {}
        self._nb_types: dict[tuple[str, str], tuple[float, float]] = {}
        self._exclusions: None | ExclusionHelper = None
        self._list = None

    def get_exclusion_helper(self):
        """Called by S* before anything else is done with this class."""
        self._exclusions = ExclusionHelper(self, self._sysstar)
        return self._exclusions

    def add(self, members, params):
        """Atom types are listed in S*, not here"""
        raise NotImplementedError

    def get_members(self, i):
        """Non bonded force, so invalid"""
        raise NotImplementedError

    def remove(self, i):
        """Non bonded force, so invalid"""
        raise NotImplementedError

    def update_params(self, i, atom_type, charge, sc_lam, sc_alpha, charge_changed):
        atom_type_id = self.use_atom_type(atom_type)
        self._force_obj.setParticleParameters(i, [atom_type_id, charge, sc_lam, sc_alpha])
        # TODO LJ type change event?
        self._sysstar.pairs._rebuild = True
        self._sysstar.cmap._rebuild = True
        if charge_changed:
            self._exclusions._rebuild = True

    def _build(self):
        self._force_obj = mm.CustomNonbondedForce(
            "step(rcut-r)*(LJ - corr + ES);"
            "LJ = (1 - sc_lambda1) * (C12(type1, type2) / rA^12 - C6(type1, type2) / rA^6) + sc_lambda1 * (C12(type1, type2) / rB^12 - C6(type1, type2) / rB^6);"
            "rA = (sc_alpha1 * sc_sigma(type1, type2)^6 * sc_lambda1^1 + r^6)^(1/6);"
            "rB = (sc_alpha1 * sc_sigma(type1, type2)^6 * (1 - sc_lambda1)^1 + r^6)^(1/6);"
            "corr = (C12(type1, type2) / rcut^12 - C6(type1, type2) / rcut^6);"
            "ES = f/epsilon_r*q1*q2 * (1/r + krf * r^2 - crf);"
            "crf = 1 / rcut + krf * rcut^2;"
            "krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {self.epsilon_r};"
            "f = 138.935458;"
            f"rcut={self.cutoff_nm};"
        )
        self._force_obj.addPerParticleParameter("type")
        self._force_obj.addPerParticleParameter("q")
        self._force_obj.addPerParticleParameter("sc_lambda")
        self._force_obj.addPerParticleParameter("sc_alpha")
        self._force_obj.setNonbondedMethod(
            mm.CustomNonbondedForce.CutoffPeriodic
        )
        self._force_obj.setCutoffDistance(self.cutoff_nm)

        for i in range(self._sysstar.len_atoms()):
            type, charge, _, sc_lam, sc_alpha = self._sysstar.get_atom_details(i)
            atom_type_id = self.use_atom_type(type)
            self._force_obj.addParticle([atom_type_id, charge, sc_lam, sc_alpha])

        for (i, j) in filter(None, self._exclusions._list):
            self._force_obj.addExclusion(i, j)

        # add LJ parameters to the system
        C6 = []
        C12 = []
        sc_sigma = []
        # i,j => type index; t1,t2 => type names
        n = len(self._used_atom_types)
        for t1, i in self._used_atom_types.items():
            for t2, j in self._used_atom_types.items():
                sigma, epsilon = (
                    self._nb_types.get((t1, t2))
                    or self._nb_types.get((t2, t1))
                    or (0., 0.)
                )
                c6 = 4 * epsilon * (sigma ** 6)
                c12 = 4 * epsilon * (sigma ** 12)
                C6.append(c6)
                C12.append(c12)
                if c6 == 0 or c12 == 0:
                    sc_sigma_value = 0.3
                else:
                    sc_sigma_value = (c12 / c6) ** (1/6)
                sc_sigma.append(sc_sigma_value)
        self._force_obj.addTabulatedFunction(
            "C6", mm.Discrete2DFunction(n, n, C6)
        )
        self._force_obj.addTabulatedFunction(
            "C12", mm.Discrete2DFunction(n, n, C12)
        )
        self._force_obj.addTabulatedFunction(
            "sc_sigma", mm.Discrete2DFunction(n, n, sc_sigma)
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

    def use_atom_type(self, atom_type):
        """Add new atom types for LJ.

        It is not a problem to call this with atom types already present.

        Returns the atom type index
        """
        id = self._used_atom_types.get(atom_type)
        if id is None:
            id = len(self._used_atom_types)
            self._used_atom_types[atom_type] = id
            self._rebuild = True
        return id

    def _additional_save(self, f):
        f.dump(self.epsilon_r)
        f.dump(self.cutoff_nm)
        f.dump(self._nb_types)
        f.dump(self._used_atom_types)

    def _additional_load(self, f):
        self.epsilon_r = f.load()
        self.cutoff_nm = f.load()
        self._nb_types = f.load()
        self._used_atom_types = f.load()


class ExclusionHelper(Force):
    """
    Just a wrapper around handling exclusions using Interaction
    handled by NonBonded

    Also handles the electrostatic self correction force

    Different to other Force objects because it is partly managed by NonBonded.

    See:
    https://manual.gromacs.org/documentation/current/reference-manual/functions/nonbonded-interactions.html
    """
    _list: list[tuple[int, int]]
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
        _, q1, _, _, _ = self._sysstar.get_atom_details(i)
        _, q2, _, _, _ = self._sysstar.get_atom_details(j)
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
            f"epsilon_r = {self._nb.epsilon_r};"
            f"f = 138.935458;"
            f"rcut={self._nb.cutoff_nm};"
        )
        self._force_obj.addPerBondParameter("q_product")
        self._force_obj.setUsesPeriodicBoundaryConditions(True)
        for i in range(self._sysstar.len_atoms()):
            _, charge, _, _, _ = self._sysstar.get_atom_details(i)
            if charge != 0:
                # self term in reaction field correction
                self.es_self_correction_add(i, i)
        for (i, j) in filter(None, self._list):
            self.es_self_correction_add(i, j)

    def add(self, i, j):
        self._list.append((i, j))
        if not self._rebuild:
            self.es_self_correction_add(i, j)
            self._sysstar._reinitialize = True
        if not self._nb._rebuild:
            self._nb._force_obj.addExclusion(i, j)
            self._sysstar._reinitialize = True

        return Interaction(self, len(self._list) - 1)

    def get_members(self, i):
        return self._list[i]

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
