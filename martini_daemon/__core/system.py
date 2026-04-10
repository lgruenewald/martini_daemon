from collections import OrderedDict
from typing import Any, Collection, Iterable, Type

import openmm as mm

from ..__rust import BondGraph, PeriodicBox
from .bonded_force import BondedForce
from .force import Force
from .molecule_type import MoleculeType


class System:
    __available_forces: dict[str, Type[BondedForce]] = {}

    def __init__(self, options=None):
        self.__names: list[str] = []
        self.__res_ids: list[int] = []
        self.__res_names: list[str] = []
        self.__types: list[str] = []
        self.__charges: list[float] = []
        self.__masses: list[float] = []
        # sc_lambda, sc_alpha
        self.__softcore: list[tuple[float, float]] = []

        # default charge, mass
        self.__atom_types: OrderedDict[str, tuple[float, float]] = OrderedDict()

        self.__context = None
        self.__system: mm.System = mm.System()
        self.__forces: dict[str, Force] = {}

        self.__interactions_by_atom: list[list[tuple[str, int]]] = []

        self.__last_resid = 0
        self.__last_force_group = 1

        self.__harmonic_constraints = None

        self.molecule_types: dict[str, MoleculeType] = {}
        self.initial_molecules: list[tuple[str, int]] = []
        self.additional_data: dict[str, Any] = options or {}

    def __assert_no_context(self):
        if self.__context is not None:
            raise ValueError(
                "This operation must be done before Context is initialized."
            )

    # ==== Initial molecules ====

    def build_initial_molecules(self):
        for name, n in self.initial_molecules:
            mol = self.molecule_types[name]
            for i in range(n):
                atoms = mol.add_atoms_to_system(self)
                mol.instantiate(self, atoms)

    # ==== ATOM METADATA ====

    def add_atom(
        self,
        name: str,
        res_name: str,
        atom_type: str,
        charge: float | None,
        mass: float | None,
        sc: tuple[float, float] = (1.0, 0.5),
    ) -> int:
        """
        Add an atom to the system.

        :returns: atom_id
        """
        self.__assert_no_context()
        defaults = self.__atom_types.get(atom_type)
        if defaults is None:
            raise ValueError(f"Unknown atom type: {atom_type}")
        if charge is None:
            charge = defaults[0]
        if mass is None:
            mass = defaults[1]
        self.__names.append(name)
        self.__res_ids.append(self.__last_resid)
        self.__res_names.append(res_name)
        self.__types.append(atom_type)
        self.__charges.append(charge)
        self.__masses.append(mass)
        self.__system.addParticle(mass)
        self.__softcore.append(sc)
        self.__interactions_by_atom.append([])
        assert (
            len(self.__names)
            == len(self.__res_ids)
            == len(self.__res_names)
            == len(self.__types)
            == len(self.__charges)
            == len(self.__masses)
            == len(self.__softcore)
        )
        self.flag_atom_add()
        return len(self.__names) - 1

    def new_residue(self):
        self.__last_resid += 1

    def atom_count(self) -> int:
        return len(self.__names)

    def get_name(self, atom_id: int) -> str:
        return self.__names[atom_id]

    def get_atom_names(self) -> list[str]:
        return self.__names

    def get_res_id(self, atom_id: int) -> int:
        return self.__res_ids[atom_id]

    def get_res_ids(self) -> list[int]:
        return self.__res_ids

    def get_res_names(self) -> list[str]:
        return self.__res_names

    def get_res_name(self, atom_id: int) -> str:
        return self.__res_names[atom_id]

    def rename(self, atom_id: int, new_name: str) -> None:
        self.__names[atom_id] = new_name

    def get_type(self, atom_id: int) -> str:
        return self.__types[atom_id]

    def get_types(self) -> list[str]:
        return self.__types

    def retype(self, atom_id: int, new_type: str) -> None:
        self.__types[atom_id] = new_type
        self.flag_atom_change(atom_id)

    def get_charge(self, atom_id: int) -> float:
        return self.__charges[atom_id]

    def get_charges(self) -> list[float]:
        return self.__charges

    def recharge(self, atom_id: int, new_charge: float) -> None:
        self.__charges[atom_id] = new_charge
        self.flag_atom_change(atom_id, True)

    def get_mass(self, atom_id: int) -> float:
        return self.__masses[atom_id]

    def get_masses(self) -> list[float]:
        return self.__masses

    def remass(self, atom_id: int, new_mass: float) -> None:
        self.__masses[atom_id] = new_mass
        self.__system.setParticleMass(atom_id, new_mass)
        self.flag_reinitialize()

    def get_sc(self, atom_id: int) -> tuple[float, float]:
        return self.__softcore[atom_id]

    def get_sc_lam(self, atom_id) -> float:
        return self.get_sc(atom_id)[0]

    def get_sc_alpha(self, atom_id) -> float:
        return self.get_sc(atom_id)[1]

    def update_sc(self, atom_id: int, sc_lam: float, sc_alpha: float) -> None:
        self.__softcore[atom_id] = (sc_lam, sc_alpha)
        self.flag_atom_change(atom_id)

    # ==== ATOM TYPE HANDLING ====

    def add_atom_type(self, atom_type: str, charge: float, mass: float) -> None:
        self.__assert_no_context()
        self.__atom_types[atom_type] = (charge, mass)

    def iterate_atom_types(self) -> Iterable[tuple[str, tuple[float, float]]]:
        return self.__atom_types.items()

    def get_atom_type(self, atom_type: str) -> None | tuple[float, float]:
        """
        Get default charge and mass for atom type.
        """
        return self.__atom_types.get(atom_type)

    # ==== STATE SYNCHRONIZATION ====

    def flag_reinitialize(self) -> None:
        """
        Mark the context associated with the system to be reinitialized.
        Doesn't do anything if the context is not yet initialized.
        """
        if self.__context is None:
            return
        self.__context.reinitialize = True

    def flag_atom_change(self, atom_id: int, change_charge: bool = False) -> None:
        for f in self.__forces.values():
            f.flag_atom_change(atom_id, change_charge)

    def flag_atom_add(self) -> None:
        for f in self.__forces.values():
            f.flag_atom_add()

    # ==== API old_helpers ====

    # PROTECTED API for __core and Forces
    def _add_mm_force(self, force: mm.Force) -> None:
        """
        Add an OpenMM force to the system. The name provided must be unique.

        Should only be called by Force/BondedForce. Should be only called if Force is in __forces.
        """
        i = self.__system.addForce(force)
        assert self.__system.getNumForces() - 1 == i
        self.flag_reinitialize()

    def _rebuild(self) -> None:
        """
        Call this before initializing or reinitializing the context.

        Should only be called by Context
        """
        for force in self.__forces.values():
            force.build()

    def _remove_mm_force(self, force: mm.Force) -> None:
        """
        Remove an OpenMM force from the system, using its unique name.

        Should only be called by Force/BondedForce.
        """
        name = force.getName()
        assert name is not None
        removed = False
        for i in range(self.__system.getNumForces()):
            if name == self.__system.getForce(i).getName():
                self.__system.removeForce(i)
                removed = True
                break
        assert removed
        self.flag_reinitialize()

    def _add_constraint(self, i, j, length) -> None:
        """
        Called by constraint force sometimes.
        """
        self.__system.addConstraint(i, j, length)

    def _add_vsite(self, atom_id: int, vsite) -> None:
        """
        Called by vsite.py.
        """
        self.__system.setVirtualSite(atom_id, vsite)

    def __del_all_constraints(self) -> None:
        """
        Used for constraints <=> harmonic bond replace.
        """
        for i in range(self.__system.getNumConstraints()-1, -1, -1):
            self.__system.removeConstraint(i)
        assert self.__system.getNumConstraints() == 0

    def toggle_constraints_as_harmonic_bonds(self, harmonic: bool) -> None:
        constraints = self.get_force("constraint")
        if constraints is None:
            return
        if harmonic:
            self.__harmonic_constraints = mm.HarmonicBondForce()
            self.__harmonic_constraints.setName("harmonic_replacement_for_constraints")
            self.__del_all_constraints()
            for _, (members, params) in constraints.iterate_bonds():
                self.__harmonic_constraints.addBond(*members, *params, 10000.)
            # this sets reinitialize to True
            self._add_mm_force(self.__harmonic_constraints)
        else:
            # we assume that only a small minimization has taken place
            # and thus no constraints were broken across pbc
            # this also flags reinitialize btw
            self._remove_mm_force(self.__harmonic_constraints)
            self.__harmonic_constraints = None
            constraints.build(must=True)

    def _set_default_pbc(self, box: PeriodicBox) -> None:
        self.__system.setDefaultPeriodicBoxVectors(box.a, box.b, box.c)

    def _bind_context(self, context):
        """
        Should only be called by context.
        """
        self.__context = context

    # PUBLIC API for users
    def add_force(self, force: Force) -> None:
        """
        Add a wrapped Force to the system.
        """
        self.__forces[force.get_name()] = force
        self.flag_reinitialize()

    @classmethod
    def provide_force(cls, force_class: Type[BondedForce]):
        """
        Register a class to the list of available forces.

        If attempting to add a new interaction, it will get instantiated and added to Forces.
        """
        assert force_class.get_name() not in cls.__available_forces.keys()
        cls.__available_forces[force_class.get_name()] = force_class

    def get_forces(self) -> Iterable[Force]:
        """
        Get an Iterator over all Forces added to System.
        """
        return self.__forces.values()

    def get_force(self, force_name: str) -> Force | None:
        """
        Get a Force by name, if present in System. None otherwise.
        """
        return self.__forces.get(force_name)

    def get_openmm_system(self):
        """
        Get the OpenMM system. Useful if you want to use Martini Daemon only as a .top file parser, and get an OpenMM
        system that you can use. Calling this will also make sure that all OpenMM Force objects get constructed
        into an appropriate state for context initialization.
        """
        self._rebuild()
        return self.__system

    def add_interaction(
        self, name: str, members: list[int], params: list[float]
    ) -> None:
        """
        Main API for getting bonds
        """
        f = self.__forces.get(name)
        if f is None:
            if name in self.__available_forces.keys():
                self.add_force(self.__available_forces[name](self))
                f = self.__forces[name]
            else:
                raise ValueError(
                    f"Internal error: {name} is not an interaction that is in System."
                )
        bond_id = f._add_bond(members, params)
        for member in members:
            assert member >= 0, f"Internal error: {member} is negative in {name}, {members}, {params}."
            self.__interactions_by_atom[member].append((name, bond_id))

    def remove_interaction(self, force_name: str, bond_id: int) -> None:
        f = self.__forces[force_name]
        members = f.get_members(bond_id)
        for member in members:
            self.__interactions_by_atom[member] = [
                (name, index)
                for name, index in self.__interactions_by_atom[member]
                if name != force_name or index != bond_id
            ]
        f._remove_bond(bond_id)

    def get_members(self, force_name: str, bond_id: int) -> list[int]:
        f = self.__forces[force_name]
        return f.get_members(bond_id)

    def get_filters(self):
        filters = set()
        for f in self.__forces.values():
            if issubclass(type(f), BondedForce):
                filters.add(f.get_name())
                for filter_ in f.filters:
                    filters.add(filter_)
        for f in self.__available_forces.values():
            if issubclass(f, BondedForce):
                filters.add(f.get_name())
                for filter_ in f.filters:
                    filters.add(filter_)
        return filters

    def break_group(self, break_group: set[int]) -> None:
        """
        Break all interactions that include all atoms in break_group.
        """
        choice = break_group.pop()
        break_group.add(choice)
        for force_name, bond_id in self.__interactions_by_atom[choice].copy():
            f = self.__forces.get(force_name)
            members = f.get_members(bond_id)
            if all(group_element in members for group_element in break_group):
                self.remove_interaction(force_name, bond_id)

    def update_group(self, update_group: set[int]) -> None:
        """
        Break all interactions that are fully contained within update_group.
        """
        for choice in update_group:
            for force_name, bond_id in self.__interactions_by_atom[choice].copy():
                f = self.__forces.get(force_name)
                members = f.get_members(bond_id)
                if all(member in update_group for member in members):
                    self.remove_interaction(force_name, bond_id)

    def get_interactions_for_atom(self, atom_id: int) -> list[tuple[str, int]]:
        """
        Get a list of force names and bond_ids that atom_id is in.
        """
        return self.__interactions_by_atom[atom_id]

    def get_number_of_degrees_of_freedom(self) -> int:
        return self.atom_count() - sum(
            f.delta_degrees_of_freedom() for f in self.__forces.values()
        )

    def collect_bonds(self, filters: list[str]) -> BondGraph:
        """
        Returns a list of bonds that matches any of the filters.

        Virtual sites will be added as each constructing particle being bonded to the virtual site.

        Multi-members interactions otherwise will be added as each member being bonded to the next.
        """
        n = self.atom_count()
        bonds = BondGraph(n)
        for force in self.get_forces():
            if not issubclass(type(force), BondedForce):
                continue
            # do we include this force
            do_force = False
            if filters is None:
                do_force = True
            else:
                for filt in filters:
                    if force.passes_filter(filt):
                        do_force = True
                        break
            if not do_force:
                continue
            for _, (members, _) in force.iterate_bonds():
                if force.passes_filter("vsite"):
                    # hardcoded special case, modeled as vsite bonded to all constructing particles
                    i = members[0]
                    for j in members[1:]:
                        if i == j:
                            continue
                        bonds.add_bond(i, j)
                else:
                    # modeled as each particle bonded to the next one
                    for i, j in zip(members[:-1], members[1:]):
                        if i == j:
                            continue
                        bonds.add_bond(i, j)
        return bonds

    def collect_bonds_for_whole(self) -> BondGraph:
        """
        Returns a bond graph that need to be made whole across the PBC before simulation
        can start. This usually includes constraints and virtual sites.
        Uses the `uses_pbc()` function of BondedForce to determine which one it is.
        """
        n = self.atom_count()
        bonds = BondGraph(n)
        for force in self.get_forces():
            # do we include this force
            if not issubclass(type(force), BondedForce) or force.uses_pbc():
                continue
            for _, (members, _) in force.iterate_bonds():
                i = members[0]
                for j in members[1:]:
                    if i == j:
                        continue
                    bonds.add_bond(i, j)
        return bonds

    def populate_neighbors(self, atoms: Collection[int], recursive=False) -> set[int]:
        """
        For a set of atoms, return a set that also contains their neighbors.

        Neighbor = shared interaction (e.g. bond, angle, exclusion...).

        If recursive, will traverse interactions recursively to get the whole molecule.
        """
        res = set()
        stack = list(atoms)
        while len(stack) > 0:
            atom = stack[-1]
            del stack[-1]
            if recursive and atom in res:
                continue
            res.add(atom)
            for force, bond_id in self.__interactions_by_atom[atom]:
                for m in self.get_force(force).get_members(bond_id):
                    assert m >= 0, f"Internal error: m is {m}, for force {force}, bond_id {bond_id}."
                    if recursive:
                        stack.append(m)
                    else:
                        res.add(m)
        return res


def register_available_force(cls):
    System.provide_force(cls)
    return cls
