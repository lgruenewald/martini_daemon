from __future__ import annotations
import openmm as mm
from openmm.unit import nanometer
from .force import Force, Interaction
from collections import OrderedDict


# mie math helpers
def mie_potential(x, m, n, sigma, epsilon):
    C = (n) / (n-m) * ((n/m) ** (m/(n-m)))
    return C * epsilon * ((sigma / x) ** n - (sigma / x) ** m)


def mie_potential_minimum(m, n, sigma):
    return ((m * sigma ** m) / (n * sigma ** n)) ** (1/(m-n))


def mie_force(x, m, n, sigma, epsilon):
    # Technically a misnomer, this is just V'(x),
    # a force with the right direction too would be -V'(x).
    C = (n) / (n-m) * ((n/m) ** (m/(n-m)))
    return C * epsilon * (
        -n * (sigma ** n) / (x ** (n+1)) +
        m * (sigma ** m) / (x ** (m+1))
    )


def shift(x, b, c, d, s):
    x = x - s
    return b * (x ** 2) / 2 + c * (x ** 3) / 6 + d * (x ** 4) / 24


def mie_shifted(x, m, n, b, c, d, s, sigma, epsilon):
    return mie_potential(x, m, n, sigma, epsilon) + shift(x, b, c, d, s)


class NonBonded(Force):
    """
    S* modular_forces implementation of non bonded forces.
    Different to other Force objects because of the interplay between the
    NB force, self exclusion force, exclusions and atomtypes.
    """

    """Non bonded parameters are for building the C6/C12 table
    values are (type1: string, type2: string), (sigma: float, epsilon: float)
    """
    _nb_types: dict[tuple[str, str], tuple[float, float]]
    _used_atom_types: OrderedDict[str, int]
    _exclusions: ExclusionHelper
    _rebuild: bool
    _force_obj: mm.Force

    def __init__(self, sysstar, nonbonded_type):
        self._sysstar = sysstar
        self._rebuild = True
        self._force_obj = None
        self._used_atom_types = {}
        self._nb_types = {}
        self._exclusions = ExclusionHelper(self, sysstar)
        self._list = None
        self.nonbonded_type = nonbonded_type

    def add(self, members, params):
        """Atom types are listed in S*, not here"""
        raise NotImplementedError

    def get_members(self, i):
        """Non bonded force, so invalid"""
        raise NotImplementedError

    def remove(self, i):
        """Non bonded force, so invalid"""
        raise NotImplementedError

    def update_params(self, i, atom_type, charge, charge_changed):
        atom_type_id = self.use_atom_type(atom_type)
        self._force_obj.setParticleParameters(i, [atom_type_id, charge])
        self._sysstar.pairs._rebuild = True
        if charge_changed:
            self._es_force._rebuild = True

    def _build(self):
        cutoff = self._sysstar.nonbonded_cutoff.value_in_unit(nanometer)
        epsilon_r = self._sysstar.epsilon_r
        if self.nonbonded_type == "default":
            self._force_obj = mm.CustomNonbondedForce(
                "step(rcut-r)*(LJ - corr + ES);"
                "LJ = (C12(type1, type2) / r^12 - C6(type1, type2) / r^6);"
                "corr = (C12(type1, type2) / rcut^12 - C6(type1, type2) / rcut^6);"
                "ES = f/epsilon_r*q1*q2 * (1/r + krf * r^2 - crf);"
                "crf = 1 / rcut + krf * rcut^2;"
                "krf = 1 / (2 * rcut^3);"
                f"epsilon_r = {epsilon_r};"
                "f = 138.935458;"
                f"rcut={cutoff};"
            )
        elif self.nonbonded_type[:3] == "mie":
            _, mie_n, mie_m = self.nonbonded_type.split("-")
            mie_n = int(mie_n)
            mie_m = int(mie_m)
            C = (mie_n) / (mie_n-mie_m) * ((mie_n/mie_m) ** (mie_m/(mie_n-mie_m)))
            self._force_obj = mm.CustomNonbondedForce(
                "step(rcut-r)*(F + S + ES);"
                f"F={C}*epsilon(type1,type2)*(term^{mie_n}-term^{mie_m});"
                "term=sigma(type1,type2)/r;"
                "S=step(switch(type1,type2)-r)*r*r*(b(type1,type2)/2+r*c(type1,type2)/6+r*r*d(type1,type2)/24);"
                "ES = f/epsilon_r*q1*q2 * (1/r + krf * r^2 - crf);"
                "crf = 1 / rcut + krf * rcut^2;"
                "krf = 1 / (2 * rcut^3);"
                f"epsilon_r = {epsilon_r};"
                "f = 138.935458;"
                f"rcut={cutoff};"
            )

        else:
            raise ValueError(f"Unknown nonbonded type: {self.nonbonded_type}")
        self._force_obj.addPerParticleParameter("type")
        self._force_obj.addPerParticleParameter("q")
        self._force_obj.setNonbondedMethod(
            mm.CustomNonbondedForce.CutoffPeriodic
        )
        self._force_obj.setCutoffDistance(
            self._sysstar.nonbonded_cutoff.value_in_unit(nanometer)
        )

        for i in range(self._sysstar.len_atoms()):
            type, charge, _ = self._sysstar.get_atom_details(i)
            atom_type_id = self.use_atom_type(type)
            self._force_obj.addParticle([atom_type_id, charge])

        for (i, j) in filter(None, self._exclusions._list):
            self._force_obj.addExclusion(i, j)

        # NOTE:
        # top_parser does sigma, epsilon -> C6, C12 parsing if non bonded type
        # is default, otherwise the values are actuall sigma and epsilon

        # add LJ parameters to the system
        C6 = []
        C12 = []
        # i,j => type index; t1,t2 => type names
        n = len(self._used_atom_types)
        for t1, i in self._used_atom_types.items():
            for t2, j in self._used_atom_types.items():
                c6, c12 = self._nb_types.get((t1, t2)) or\
                            self._nb_types.get((t2, t1)) or (0., 0.)
                C6.append(c6)
                C12.append(c12)
        if self.nonbonded_type == "default":
            self._force_obj.addTabulatedFunction(
                "C6", mm.Discrete2DFunction(n, n, C6)
            )
            self._force_obj.addTabulatedFunction(
                "C12", mm.Discrete2DFunction(n, n, C12)
            )
        elif self.nonbonded_type[:3] == "mie":
            # mie_n, mie_m defined above conditionally with the same condition
            # adjusted sigmas
            sigma = []
            # switch+shift param quadratic
            b = []
            # switch+shift param cubic
            c = []
            # switch+shift param quartic
            d = []
            # switch point
            switch = []
            for i in range(len(C6)):
                r_min_n_m = mie_potential_minimum(mie_m, mie_n, C6[i])
                r_min = mie_potential_minimum(6, 12, C6[i])
                sigma.append(C6[i] * r_min / r_min_n_m)

                # switching constant calculation
                F_cutoff = mie_potential(cutoff, mie_m, mie_n, sigma[i], C12[i])
                F_prime_cutoff = mie_force(cutoff, mie_m, mie_n, sigma[i], C12[i])
                diff = cutoff - r_min
                # quadratic-cubic shift+switch:
#                b.append(2 * F_prime_cutoff / diff - 6 * F_cutoff / (diff ** 2))
#                c.append(12 * F_cutoff / (diff ** 3) - 6 * F_prime_cutoff / (diff ** 2))
                # cubic-quartic shift+switch:
                b.append(0.)
                c.append(6 * F_prime_cutoff / diff ** 2 - 24 * F_cutoff / diff ** 3)
                d.append(72 / diff ** 4 * F_cutoff - 24 * F_prime_cutoff / diff ** 3)
                switch.append(r_min)
                # quadratic-cubic shift only
                # b.append(2 * F_prime_cutoff / cutoff - 6 * F_cutoff / (cutoff ** 2))
                # c.append(12 * F_cutoff / (cutoff ** 3) - 6 * F_prime_cutoff / (cutoff ** 2))
                # d.append(0.)
                # switch.append(0.)
            assert all([n * n == len(x) for x in [sigma, C12, b, c, d, switch]])
            self._force_obj.addTabulatedFunction(
                "sigma", mm.Discrete2DFunction(n, n, sigma)
            )
            self._force_obj.addTabulatedFunction(
                "epsilon", mm.Discrete2DFunction(n, n, C12)
            )
            self._force_obj.addTabulatedFunction(
                "b", mm.Discrete2DFunction(n, n, b)
            )
            self._force_obj.addTabulatedFunction(
                "c", mm.Discrete2DFunction(n, n, c)
            )
            self._force_obj.addTabulatedFunction(
                "d", mm.Discrete2DFunction(n, n, d)
            )
            self._force_obj.addTabulatedFunction(
                "switch", mm.Discrete2DFunction(n, n, switch)
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

    def get_exclusion_helper(self):
        return self._exclusions

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

    def is_instance(self, filter):
        return False

    def _additional_save(self, f):
        f.dump(self._nb_types)
        f.dump(self._used_atom_types)

    def _additional_load(self, f):
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
    # TODO is a dict based approach better, then duplicate additions of
    # exclusions could return the same one
    # does returning an existing exclusion create any problems? what if that
    # exclusion gets removed in a reaction?
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
        _, q1, _ = self._sysstar.get_atom_details(i)
        _, q2, _ = self._sysstar.get_atom_details(j)
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
        for i in range(self._sysstar.len_atoms()):
            _, charge, _ = self._sysstar.get_atom_details(i)
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

    def is_instance(self, filter):
        return False
