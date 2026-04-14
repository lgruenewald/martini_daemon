from __future__ import annotations

from collections import OrderedDict

import openmm as mm

from ..__core import BondedForce, Force


class NonBonded(Force):
    """A special force. It is special, because it is passed as an argument to Simulation rather than constructed
    based on .top files. It also is responsible for generating the Exclusion Helper. Simulation handles adding it
    to System, which does not treat it in a special way. The Exclusion Helper's name "exclusion" is also special,
    since molecule types hardcode it, since exclusion handling is its own whole thing.

    Exclusions are kept up to date in this manner:
    - when build()-ing, read all exclusions from exclusion helper.
    - when already built, exclusion helper calls add_exclusion(), which will add it, then flag reinit.
    - when removing exclusions, just remove the whole force and rebuild, we don't keep track of indices and it's not a bottleneck.
    """

    @classmethod
    def is_coupling(cls):
        return False

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "nonbonded"

    def __init__(self, system):
        super().__init__(system)
        self.epsilon_r = system.additional_data.get("epsilon_r")
        self.cutoff_nm = system.additional_data.get("cutoff")
        self.__atom_types: OrderedDict[str, int] = OrderedDict()
        self.__exclusions: None | ExclusionHelper = None

    def get_exclusion_helper(self):
        """Called by Simulation when adding the NonBonded force to the system.

        ExclusionHelper is responsible for:
        - adding and removing exclusions from NonBonded and flagging reinitialize when that happens.
        """
        self.__exclusions = ExclusionHelper(self.system, self)
        return self.__exclusions

    def flag_atom_change(self, atom_id: int, change_charge: bool = False) -> None:
        if self.force is None:
            return
        type_ = self.__atom_types[self.system.get_type(atom_id)]
        charge = self.system.get_charge(atom_id)
        sc_lam, sc_alpha = self.system.get_sc(atom_id)
        assert isinstance(self.force, mm.CustomNonbondedForce)
        self.force.setParticleParameters(atom_id, [type_, charge, sc_lam, sc_alpha])
        self.system.flag_reinitialize()

    def flag_atom_add(self):
        self._destroy()

    def add_exclusion(self, i: int, j: int) -> None:
        if self.force is not None:
            assert isinstance(self.force, mm.CustomNonbondedForce)
            self.force.addExclusion(i, j)
            self.system.flag_reinitialize()

    def flag_remove_exclusion(self) -> None:
        self._destroy()

    def _set_force_obj(self):
        if self.__exclusions is None:
            raise AssertionError("Must call get_exclusion_helper() first!")
        assert isinstance(self.__exclusions, ExclusionHelper)
        self.force = mm.CustomNonbondedForce(
            "(LJ - corr + ES);"
            "LJ = (1 - sc_lambda1) * (C12 / rA^2 - C6 / rA) + sc_lambda1 * (C12 / rB^2 - C6 / rB);"
            "rA = (sc_alpha1 * sigma(type1, type2)^6 * sc_lambda1 + r^6);"
            "rB = (sc_alpha1 * sigma(type1, type2)^6 * (1 - sc_lambda1) + r^6);"
            "corr = (C12 / rcut^12 - C6 / rcut^6);"
            "C6 = 4 * epsilon(type1, type2) * sigma(type1, type2)^6;"
            "C12 = 4 * epsilon(type1, type2) * sigma(type1, type2)^12;"
            "ES = f/epsilon_r*q1*q2 * (1/r + krf * r^2 - crf);"
            "crf = 1 / rcut + krf * rcut^2;"
            "krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {self.epsilon_r};"
            "f = 138.935458;"
            f"rcut = {self.cutoff_nm};"
        )
        self.force.addPerParticleParameter("type")
        self.force.addPerParticleParameter("q")
        self.force.addPerParticleParameter("sc_lambda")
        self.force.addPerParticleParameter("sc_alpha")
        self.force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        self.force.setCutoffDistance(self.cutoff_nm)

        for i, (type_name, _) in enumerate(self.system.iterate_atom_types()):
            if self.__atom_types.get(type_name) is not None:
                # this shouldn't happen
                assert i == self.__atom_types[type_name], (
                    "Internal error: atom types changed"
                )
            self.__atom_types[type_name] = i

        for atom_id in range(self.system.num_atoms()):
            type_ = self.__atom_types[self.system.get_type(atom_id)]
            charge = self.system.get_charge(atom_id)
            sc_lam, sc_alpha = self.system.get_sc(atom_id)
            self.force.addParticle([type_, charge, sc_lam, sc_alpha])

        for _, (members, _) in self.__exclusions.iterate_bonds():
            self.force.addExclusion(*members)

        # add LJ parameters to the system
        sigmas = []
        epsilons = []
        # i,j => type index; t1,t2 => type names
        n = len(self.__atom_types)
        nb_types: dict[tuple[str, str], tuple[float, float]] = (
            self.system.additional_data.get("nb_types")
        )
        for t1, i in self.__atom_types.items():
            for t2, j in self.__atom_types.items():
                sigma, epsilon = (
                    nb_types.get((t1, t2)) or nb_types.get((t2, t1)) or (0.0, 0.0)
                )
                sigmas.append(sigma)
                epsilons.append(epsilon)
        self.force.addTabulatedFunction("sigma", mm.Discrete2DFunction(n, n, sigmas))
        self.force.addTabulatedFunction(
            "epsilon", mm.Discrete2DFunction(n, n, epsilons)
        )


class ExclusionHelper(BondedForce):
    """The force behind the force name "exclusion". Instantiated by NonBonded.get_exclusion_helper().
    Added to the system by Simulation. Each NonBonded force implementation should provide its own.
    Also handles the electrostatic self correction force

    Different to other BondedForce objects because it is partly managed by NonBonded / partly manages NonBonded.

    For electrostatic self correction see:
    https://manual.gromacs.org/documentation/current/reference-manual/functions/nonbonded-interactions.html
    """

    def _add_to_force(
        self, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        assert isinstance(force, mm.CustomBondForce)
        self.es_self_correction_add(force, *members)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @staticmethod
    def uses_pbc() -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    def should_build(self) -> bool:
        # should build also when there is no exclusions, so that self corrections work
        return self.force is None

    def _set_force_obj(self) -> None:
        self.force = mm.CustomBondForce(
            f"step(rcut-r) * ES;"
            f"ES = f*q_product/epsilon_r * (krf * r^2 - crf);"
            f"crf = 1 / rcut + krf * rcut^2;"
            f"krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {self.epsilon_r};"
            f"f = 138.935458;"
            f"rcut={self.cutoff_nm};"
        )
        self.force.addPerBondParameter("q_product")

        # trigger rebuild! if a charge changes
        for i in range(self.system.num_atoms()):
            charge = self.system.get_charge(i)
            if charge != 0:
                # self term in reaction field correction
                self.es_self_correction_add(self.force, i, i)

    @classmethod
    def get_name(cls) -> str:
        return "exclusion"

    def __init__(self, system, nb):
        super().__init__(system)
        self.__nb = nb
        self.epsilon_r = system.additional_data.get("epsilon_r")
        self.cutoff_nm = system.additional_data.get("cutoff")

    def es_self_correction_add(self, force: mm.CustomBondForce, i, j):
        q1 = self.system.get_charge(i)
        q2 = self.system.get_charge(j)
        q_prod = q1 * q2
        if i == j:
            q_prod *= 0.5
        if q_prod != 0:
            force.addBond(i, j, [q_prod])

    def _add_bond(self, members: list[int], params: list[float]) -> int:
        res = super()._add_bond(members, params)
        self.__nb.add_exclusion(*members)
        return res

    def _remove_bond(self, bond_id: int) -> None:
        super()._remove_bond(bond_id)
        self.__nb.flag_remove_exclusion()

    def flag_atom_add(self):
        self._destroy()

    def flag_atom_change(self, atom_id: int, change_charge: bool = False) -> None:
        if change_charge:
            self._destroy()
