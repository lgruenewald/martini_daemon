import openmm as mm
from typing import Iterable
from collections import OrderedDict
from .force import Force

class System:

    def __init__(self):
        self.__names: list[str] = []
        self.__res_ids: list[int] = []
        self.__res_names: list[str] = []
        self.__types: list[str] = []
        self.__charges: list[float] = []
        self.__masses: list[float] = []
        self.__softcore: list[tuple[float, float]] = []

        # default charge, mass
        self.__atom_types: OrderedDict[str, tuple[float, float]] = OrderedDict()

        self.__context = None
        self.__system: mm.System = mm.System()
        self.__forces: dict[str, Force] = {}

        self.__interactions_by_atom: list[list[tuple[str, int]]]

        self.__last_resid = 0
        self.__last_force_group = 1

    def __assert_no_context(self):
        if self.__context is not None:
            raise ValueError("This operation must be done before Context is initialized.")

    # ==== ATOM METADATA ====

    def add_atom(
        self, name: str, res_name: str, atom_type: str,
        charge: float | None, mass: float | None, sc: tuple[float, float] = (1., 0.5)
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
        assert (
            len(self.__names) == len(self.__res_ids) == len(self.__res_names) == len(self.__types)
            == len(self.__charges) == len(self.__masses) == len(self.__softcore)
        )
        return len(self.__names) - 1

    def new_residue(self):
        self.__last_resid += 1

    def atom_count(self) -> int:
        return len(self.__names)

    def get_name(self, atom_id: int) -> str:
        return self.__names[atom_id]

    def rename(self, atom_id: int, new_name: str) -> None:
        self.__names[atom_id] = new_name

    def get_type(self, atom_id: int) -> str:
        return self.__types[atom_id]

    def retype(self, atom_id: int, new_type: str) -> None:
        self.__types[atom_id] = new_type
        self.flag_nonbonded(atom_id)

    def get_charge(self, atom_id: int) -> float:
        return self.__charges[atom_id]

    def recharge(self, atom_id: int, new_charge: float) -> None:
        self.__charges[atom_id] = new_charge
        self.flag_nonbonded(atom_id, True)

    def get_mass(self, atom_id: int) -> float:
        return self.__masses[atom_id]

    def remass(self, atom_id: int, new_mass: float) -> None:
        self.__masses[atom_id] = new_mass
        self.__system.setParticleMass(atom_id, new_mass)
        self.flag_reinitialize()

    def get_sc(self, atom_id: int) -> tuple[float, float]:
        return self.__softcore[atom_id]

    def update_sc(self, atom_id: int, new_sc: tuple[float, float]) -> None:
        self.__softcore[atom_id] = new_sc
        self.flag_nonbonded(atom_id)

    # ==== ATOM TYPE HANDLING ====

    def add_atom_type(self, atom_type: str, charge: float, mass: float) -> None:
        self.__assert_no_context()
        self.__atom_types[atom_type] = (charge, mass)

    # ==== STATE SYNCHRONIZATION ====

    def flag_reinitialize(self) -> None:
        """
        Mark the context associated with the system to be reinitialized.
        Doesn't do anything if the context is not yet initialized.
        """
        if self.__context is None:
            return
        self.__context.reinitialize = True

    def flag_nonbonded(self, atom_id: int, change_charge: bool = False) -> None:
        if self.__context is None:
            return
        # 1. flag all NonBonded forces
        # nonbonded forces should flag their own exclusion helpers
        raise NotImplementedError

    # ==== API helpers ====

    def add_force(self, force: Force) -> None:
        self.__forces[force.get_name()] = force
        force.build()

    def get_forces(self) -> Iterable[Force]:
        """
        WARNING! Assume Force provided is read-only.
        """
        return self.__forces.values()

    def add_mm_force(self, force: mm.Force) -> None:
        self.__system.addForce(force)
        self.flag_reinitialize()

    def rebuild(self) -> None:
        for force in self.__forces.values():
            force.build()

    def remove_mm_force(self, force: mm.Force) -> None:
        name = force.getName()
        assert name is not None and name in self.__forces
        for i in range(self.__system.getNumForces()):
            if name == self.__system.getForce(i).getName():
                self.__system.removeForce(i)
        self.flag_reinitialize()

    def bind_context(self, context):
        self.__context = context

    def get_openmm_system(self):
        return self.__system

    def add_bond(self, name: str, members: list[int], params: list[float]) -> None:
        raise NotImplementedError

    def remove_bond(self, name: str, bond_id: int) -> None:
        raise NotImplementedError

    def get_interactions_for_atom(self, atom_id: int) -> list[tuple[str, int]]:
        raise NotImplementedError

    def get_number_of_degrees_of_freedom(self) -> int:
        raise NotImplementedError
