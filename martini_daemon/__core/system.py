import itertools
from collections import OrderedDict
from collections.abc import Iterable
from typing import Any

import openmm as mm

from ..__rust import BondGraph, PeriodicBox
from .bonded_force import BondedForce
from .force import Force
from .molecule_type import MoleculeType


class System:
    __available_forces: dict[str, type[BondedForce]] = {}

    def __init__(self, options: None | dict[str, Any] = None) -> None:
        """
        Create a Martini Daemon System.

        * Contains all the topology information.
        * Wraps OpenMM Force and System objects, and keeps them up to date based on its own topology information.
        * Provides an interface for the rest of Martini Daemon to manipulate topology information.
        * Stores global per-system information in fields, such as molecule_types, initial_molecules and additional_data.

        """
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

    def __assert_no_context(self) -> None:
        if self.__context is not None:
            raise ValueError(
                "This operation must be done before Context is initialized."
            )

    # ==== Initial molecules ====

    def build_initial_molecules(self) -> None:
        """
        Add atoms and interactions to the system based on the initial_molecules attribute.

        Note: automatically called during parsing.
        """
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
    ) -> int:
        """
        Add an atom to the system.

        Note: residue ID is assigned by system and can be bumped using System.new_residue().

        :param name: Atom name
        :param res_name: Residue name this atom belongs to
        :param atom_type: Atom type for Non-Bonded interactions
        :param charge: Atom (partial) charge
        :param mass: Atom mass, in atomic units
        :return: atom_id
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
        self.__softcore.append((1.0, 0.5))
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
        self.__flag_atom_add()
        return len(self.__names) - 1

    def new_residue(self) -> None:
        """Bump the current residue ID, used for adding new atoms to the system."""
        self.__last_resid += 1

    def num_atoms(self) -> int:
        """Return the current number of atoms in the system."""
        return len(self.__names)

    def get_name(self, atom_id: int) -> str:
        """Get the atom name for atom_id."""
        return self.__names[atom_id]

    def get_atom_names(self) -> list[str]:
        """
        Return all atom names.

        Warning! Assume returned list is read-only.
        """
        return self.__names

    def get_res_id(self, atom_id: int) -> int:
        """Get the residue number of atom_id."""
        return self.__res_ids[atom_id]

    def get_res_ids(self) -> list[int]:
        """
        Return all residue ids.

        Warning! Assume returned list is read-only.
        """
        return self.__res_ids

    def get_res_names(self) -> list[str]:
        """
        Return all residue names.

        Warning! Assume returned list is read-only.
        """
        return self.__res_names

    def get_res_name(self, atom_id: int) -> str:
        """Get the residue name of atom_id."""
        return self.__res_names[atom_id]

    def rename(self, atom_id: int, new_name: str) -> None:
        """Change the atom name of atom_id to new_name."""
        self.__names[atom_id] = new_name

    def get_type(self, atom_id: int) -> str:
        """Get the atom type of atom_id."""
        return self.__types[atom_id]

    def get_types(self) -> list[str]:
        """
        Return all atom types.

        Warning! Assume returned list is read-only.
        """
        return self.__types

    def retype(self, atom_id: int, new_type: str) -> None:
        """
        Change the atom type of atom_id to new_type.

        Informs all forces to let them automatically update it.
        """
        self.__types[atom_id] = new_type
        self.__flag_atom_change(atom_id, False)

    def get_charge(self, atom_id: int) -> float:
        """Get the (partial) charge of atom_id."""
        return self.__charges[atom_id]

    def get_charges(self) -> list[float]:
        """
        Return all (partial) charges.

        Warning! Assume returned list is read-only.
        """
        return self.__charges

    def recharge(self, atom_id: int, new_charge: float) -> None:
        """Change the (partial) charge of atom_id to new_charge."""
        self.__charges[atom_id] = new_charge
        self.__flag_atom_change(atom_id, True)

    def get_mass(self, atom_id: int) -> float:
        """Get the mass of atom_id, in atomic units."""
        return self.__masses[atom_id]

    def get_masses(self) -> list[float]:
        """
        Return all atom masses, in atomic units.

        Warning! Assume returned list is read-only.
        """
        return self.__masses

    def remass(self, atom_id: int, new_mass: float) -> None:
        """
        Change the mass of atom_id to new_mass.

        Informs all forces to let them automatically update it.
        """
        self.__masses[atom_id] = new_mass
        self.__system.setParticleMass(atom_id, new_mass)
        self.flag_reinitialize()

    def get_sc(self, atom_id: int) -> tuple[float, float]:
        """
        Get the current soft-core parameters of atom_id.

        Soft-core parameters are additional per-atom data for the Non-Bonded force, used by Local Minimization.
        """
        return self.__softcore[atom_id]

    def get_sc_lam(self, atom_id: int) -> float:
        """Get the lambda soft core parameter for atom_id."""
        return self.get_sc(atom_id)[0]

    def get_sc_alpha(self, atom_id: int) -> float:
        """Get the alpha soft core parameter for atom_id."""
        return self.get_sc(atom_id)[1]

    def update_sc(self, atom_id: int, sc_lam: float, sc_alpha: float) -> None:
        """
        Update the soft-core parameters for atom_id.

        :param atom_id: The atom id.
        :param sc_lam: The lambda soft core parameter.
        :param sc_alpha: The alpha soft core parameter.

        Note, the local minimizer will automatically do this based on soft core directives in reaction templates.
        """
        self.__softcore[atom_id] = (sc_lam, sc_alpha)
        self.__flag_atom_change(atom_id)

    # ==== ATOM TYPE HANDLING ====

    def add_atom_type(self, atom_type: str, charge: float, mass: float) -> None:
        """
        Add a new atom type to the system.

        :param atom_type: The new NonBonded atom type.
        :param charge: The default (partial) charge for this type.
        :param mass: The default mass for this type.
        """
        self.__assert_no_context()
        self.__atom_types[atom_type] = (charge, mass)

    def iterate_atom_types(self) -> Iterable[tuple[str, tuple[float, float]]]:
        """
        Get an iterator over all atom types.

        :return: An iterator over tuples of (atom_type, tuples of (charge, mass)).
        """
        return self.__atom_types.items()

    def get_atom_type(self, atom_type: str) -> None | tuple[float, float]:
        """
        Get default charge and mass for atom type.

        :param atom_type: The atom type.
        :return: The default charge and mass for atom type as a tuple.
        """
        return self.__atom_types.get(atom_type)

    # ==== STATE SYNCHRONIZATION ====

    def flag_reinitialize(self) -> None:
        """
        Mark the context associated with the system to be reinitialized.

        Doesn't do anything if the context is not yet initialized.

        Called by Forces in this system, when their contents change.
        """
        if self.__context is None:
            return
        self.__context.reinitialize = True

    def __flag_atom_change(self, atom_id: int, change_charge: bool = False) -> None:
        for f in self.__forces.values():
            f.flag_atom_change(atom_id, change_charge)

    def __flag_atom_add(self) -> None:
        for f in self.__forces.values():
            f.flag_atom_add()

    # ==== global exclusion sync ====
    def flag_remove_exclusion(self) -> None:
        """
        Tell all forces that an exclusion was removed.

        Called by Exclusion Helper.
        """
        for f in self.__forces.values():
            f.flag_remove_exclusion()

    def flag_add_exclusion(self, i: int, j: int) -> None:
        """
        Tell all forces that an exclusion was added.

        Called by Exclusion Helper.
        """
        for f in self.__forces.values():
            f.flag_add_exclusion(i, j)

    def iterate_exclusions(self) -> Iterable[tuple[int, int]]:
        excl = self.get_force("exclusion")
        assert excl is not None
        assert isinstance(excl, BondedForce)
        # FIXME: don't assume exclusion like this
        return ((m[0], m[1]) for (_, (m, _)) in excl.iterate_bonds())

    # ==== API old_helpers ====

    # PROTECTED API for __core and Forces
    def _add_mm_force(self, force: mm.Force) -> None:
        """
        Add an OpenMM force to the system. The name provided must be unique.

        Should only be called by Force/BondedForce. Should be only called if Force is already in __forces.

        User-defined forces inheriting Force/BondedForce need not worry about calling this.
        This should be considered private, changing this is not a breaking change.
        """
        i = self.__system.addForce(force)
        assert self.__system.getNumForces() - 1 == i
        self.flag_reinitialize()

    def _rebuild(self) -> None:
        """
        Call this before initializing or reinitializing the context.

        Protected because it should only be called by Context.
        This should be considered private, changing this is not a breaking change.
        """
        for force in self.__forces.values():
            force.build()

    def _remove_mm_force(self, force: mm.Force) -> None:
        """
        Remove an OpenMM force from the system, using its unique name.

        Should only be called by Force/BondedForce.

        User-defined forces inheriting Force/BondedForce need not worry about calling this.
        This should be considered private, changing this is not a breaking change.
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

    def _add_constraint(self, i: int, j: int, length: float) -> None:
        """
        Add a constraint to the OpenMM system.

        Must be called by forces implementing constraints only.

        Protected, because only Forces should call this.
        """
        self.__system.addConstraint(i, j, length)

    def _add_vsite(self, atom_id: int, vsite: mm.VirtualSite) -> None:
        """
        Add a virtual site to the OpenMM system.

        Must be called by vsite.py only.

        Protected, because only Forces/VSites should call this.
        """
        self.__system.setVirtualSite(atom_id, vsite)

    def __del_all_constraints(self) -> None:
        """
        Temporarily remove all constraints from the system.

        Used for constraints <=> harmonic bond replace.
        """
        for i in range(self.__system.getNumConstraints() - 1, -1, -1):
            self.__system.removeConstraint(i)
        assert self.__system.getNumConstraints() == 0

    def _set_default_pbc(self, box: PeriodicBox) -> None:
        """
        Set the default PBC for the OpenMM system.

        Called by context to set the default periodic box in System to it, so context can be constructed.

        Protected, but only context should call it.
        Not a part of the public API, changing this is non-breaking.
        """
        self.__system.setDefaultPeriodicBoxVectors(box.a, box.b, box.c)

    # to avoid circular dependency, no typing
    def _bind_context(self, context) -> None:  # noqa: ANN001
        """
        Set the reference for context.

        Should only be called by context.

        Protected, but only context should call it.
        Not a part of the public API, changing this is non-breaking.
        """
        self.__context = context

    def _get_context(self) -> mm.Context | None:
        """
        Get the raw OpenMM context.

        Protected, because this method currently exists to support UpdateParametersInContext
        as a potential future optimization. A better way may be provided in the future.
        """
        return self.__context._get_context() if self.__context is not None else None

    # PUBLIC API for users and custom Forces
    def toggle_constraints_as_harmonic_bonds(self, harmonic: bool) -> None:
        """
        Toggle between using harmonic bonds and constraints.

        Assumes constraints have the name "constraint".

        Note: automatically toggled on/off by LocalMinimizer.

        :param harmonic: whether to toggle on or off.
        """
        constraints = self.get_force("constraint")
        if constraints is None:
            return
        assert isinstance(constraints, BondedForce)
        if harmonic:
            self.__harmonic_constraints = mm.HarmonicBondForce()
            self.__harmonic_constraints.setName("harmonic_replacement_for_constraints")
            self.__del_all_constraints()
            for _, (members, params) in constraints.iterate_bonds():
                assert len(members) == 2
                assert len(params) == 1
                self.__harmonic_constraints.addBond(
                    members[0], members[1], params[0], 10000.0
                )
            # this sets reinitialize to True
            self._add_mm_force(self.__harmonic_constraints)
        elif self.__harmonic_constraints is not None:
            # we assume that only a small minimization has taken place
            # and thus no constraints were broken across pbc
            # this also flags reinitialize btw
            self._remove_mm_force(self.__harmonic_constraints)
            self.__harmonic_constraints = None
            constraints.build(must=True)

    def add_force(self, force: Force) -> None:
        """Add a Martini Daemon Force to the system."""
        self.__forces[force.get_name()] = force
        self.flag_reinitialize()

    @classmethod
    def provide_force(cls, force_class: type[BondedForce]) -> None:
        """
        Register a class to the list of available forces.

        If attempting to add a new interaction with its name, it will get instantiated and added to Forces.
        """
        assert force_class.get_name() not in cls.__available_forces
        cls.__available_forces[force_class.get_name()] = force_class

    def get_forces(self) -> Iterable[Force]:
        """Get an Iterator over all Forces added to System."""
        return self.__forces.values()

    def get_force(self, force_name: str) -> Force | None:
        """Get a Force by name, if present in System. None otherwise."""
        return self.__forces.get(force_name)

    def get_openmm_system(self) -> mm.System:
        """
        Get the OpenMM system.

        Useful if you want to use Martini Daemon only as a .top file parser, and get an OpenMM
        system that you can use. Calling this will also make sure that all OpenMM Force objects get constructed
        into an appropriate state for context initialization.
        """
        self._rebuild()
        return self.__system

    def add_interaction(
        self, name: str, members: list[int], params: list[float]
    ) -> None:
        """
        Add a new interaction to the System, such as bond, angle...

        Note: automatically called by MoleculeType instances or the modification algorithm as appropriate.

        :param name: name of the force.
        :param members: atom indices used for the interaction.
        :param params: raw parameters, as unwrapped from the parser for the interaction. Will call _parse() on them
            first.
        """
        f = self.__forces.get(name)
        if f is None:
            if name in self.__available_forces:
                self.add_force(self.__available_forces[name](self))
                f = self.__forces[name]
            else:
                raise ValueError(
                    f"Provided force {name} is not an interaction that is in System, nor is it available."
                )
        if not isinstance(f, BondedForce):
            raise ValueError(
                f"Provided force {name} does not correspond to a bonded force."
            )
        bond_id = f._add_bond(members, params)
        for member in members:
            assert member >= 0, (
                f"Internal error: {member} is negative in {name}, {members}, {params}."
            )
            self.__interactions_by_atom[member].append((name, bond_id))

    def remove_interaction(self, force_name: str, bond_id: int) -> None:
        """
        Remove an interaction from the System.

        :param force_name: name of the force.
        :param bond_id: interaction index.
        """
        members = self.get_members(force_name, bond_id)
        for member in members:
            self.__interactions_by_atom[member] = [
                (name, index)
                for name, index in self.__interactions_by_atom[member]
                if name != force_name or index != bond_id
            ]
        # self.get_members() already throws nice exception messages
        f = self.__forces[force_name]
        assert isinstance(f, BondedForce)
        f._remove_bond(bond_id)

    def get_members(self, force_name: str, bond_id: int) -> list[int]:
        """
        Get which members constitute an interaction.

        :param force_name: name of the force.
        :param bond_id: interaction index.
        """
        f = self.__forces.get(force_name)
        if f is None:
            raise ValueError(f"No force was found with name {force_name}.")
        if not isinstance(f, BondedForce):
            raise ValueError(f"Force with name {force_name} is not a bonded force.")
        return f.get_members(bond_id)

    def get_filters(self) -> set[str]:
        """
        Return all the possible BondedForce filters that can be used.

        Corresponds to the valid interaction filters used in the graph matching algorithm.
        """
        filters = set()
        # ALREADY INSTANTIATED
        for f in self.__forces.values():
            if isinstance(f, BondedForce):
                filters.add(f.get_name())
                for filter_ in f.filters:
                    filters.add(filter_)
        # STUFF THAT IS STILL TYPES
        for f in self.__available_forces.values():
            if issubclass(f, BondedForce):
                filters.add(f.get_name())
                for filter_ in f.filters:
                    filters.add(filter_)
        return filters

    def break_group(self, break_group: set[int]) -> None:
        """Break all interactions that include all atoms in break_group."""
        choice = break_group.pop()
        break_group.add(choice)
        for force_name, bond_id in self.__interactions_by_atom[choice].copy():
            f = self.__forces.get(force_name)
            assert isinstance(f, BondedForce)
            members = f.get_members(bond_id)
            if all(group_element in members for group_element in break_group):
                self.remove_interaction(force_name, bond_id)

    def update_group(self, update_group: set[int]) -> None:
        """Break all interactions that are fully contained within update_group."""
        for choice in update_group:
            for force_name, bond_id in self.__interactions_by_atom[choice].copy():
                f = self.__forces.get(force_name)
                assert isinstance(f, BondedForce)
                members = f.get_members(bond_id)
                if all(member in update_group for member in members):
                    self.remove_interaction(force_name, bond_id)

    def get_interactions_for_atom(self, atom_id: int) -> list[tuple[str, int]]:
        """Get a list of force names and bond_ids that atom_id is in."""
        return self.__interactions_by_atom[atom_id]

    def get_number_of_degrees_of_freedom(self) -> int:
        """
        Get the total number of degrees of freedom in the system.

        By default, 3*N, but some forces, such as virtual sites, constraints, or center of mass motion removal
        will decrease it.
        """
        return self.num_atoms() * 3 - sum(
            f.delta_degrees_of_freedom() for f in self.__forces.values()
        )

    def collect_bonds(self, filters: list[str]) -> BondGraph:
        """
        Return a list of bonds that matches any of the filters.

        Virtual sites will be added as each constructing particle being bonded to the virtual site.

        Multi-members interactions otherwise will be added as each member being bonded to the next.
        """
        n = self.num_atoms()
        bonds = BondGraph(n)
        for force in self.get_forces():
            if not isinstance(force, BondedForce):
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
                    for i, j in itertools.pairwise(members):
                        if i == j:
                            continue
                        bonds.add_bond(i, j)
        return bonds

    def collect_bonds_for_whole(self) -> BondGraph:
        """
        Return a bond graph that need to be made whole across the PBC before simulation can start.

        This usually includes constraints and virtual sites.
        Uses the `uses_pbc()` function of BondedForce to determine which one it is.
        """
        n = self.num_atoms()
        bonds = BondGraph(n)
        for force in self.get_forces():
            # do we include this force
            if not isinstance(force, BondedForce) or force.uses_pbc():
                continue
            for _, (members, _) in force.iterate_bonds():
                i = members[0]
                for j in members[1:]:
                    if i == j:
                        continue
                    bonds.add_bond(i, j)
        return bonds

    def populate_neighbors(
        self, atoms: Iterable[int], recursive: bool = False
    ) -> set[int]:
        """
        For a set of atoms, return a set that also contains their neighbors.

        Neighbors are determined based on shared interactions (e.g. bond, angle, exclusion...).
        If recursive, will traverse interactions recursively to get the whole molecule.

        :param atoms: Atom indices to determine the neighbors of.
        :param recursive: Whether to recursively get all reachable neighbors starting at atoms.
        :return: A new set of atoms, including atoms and their neighbors.
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
                f = self.get_force(force)
                assert isinstance(f, BondedForce)
                for m in f.get_members(bond_id):
                    assert m >= 0, (
                        f"Internal error: m is {m}, for force {force}, bond_id {bond_id}."
                    )
                    if recursive:
                        stack.append(m)
                    else:
                        res.add(m)
        return res


def register_available_force[T: type[BondedForce]](cls: T) -> T:
    """Add an available force to all System objects."""
    System.provide_force(cls)
    return cls
