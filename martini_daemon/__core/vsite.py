from abc import ABCMeta, abstractmethod

import openmm as mm

from .bonded_force import BondedForce
from .system import System


class VirtualSite(BondedForce, metaclass=ABCMeta):
    def __init__(self, system: System) -> None:
        """
        Override this class to create custom virtual sites.

        Base class for virtual sites, based on BondedForce.

        :param system: Martini Daemon System object.
        """
        super().__init__(system)
        self._built = False

    @abstractmethod
    def _make_vsite(
        self, vid: int, other: list[int], params: list[float]
    ) -> mm.VirtualSite:
        """
        Create an OpenMM virtual site object with given parameters.

        Overridden by child classes.

        :param vid: the atom index, which will be the virtual site, as a reference.
        :param other: list of constructing atoms.
        :param params: list of parameters.
        :return: a virtual site object.

        Note: Should not add it to system yet, should only construct the object.
        """
        raise NotImplementedError

    def build(self, must: bool = False) -> None:
        """
        Add all virtual sites to the system.

        Overrides Force's build for virtual sites, since virtual sites are not built like forces in OpenMM.
        """
        if self._built:
            return
        self._built = True
        for index, (members, params) in self.iterate_bonds():
            vid, *other = members
            self.system._add_vsite(vid, self._make_vsite(vid, other, params))

    def delta_degrees_of_freedom(self) -> int:
        """Pre-defined to return three times the number of virtual sites."""
        return self.num_bonds() * 3

    def _add_bond(self, members: list[int], params: list[float]) -> int:
        """
        Add a new virtual site.

        Same as in BondedForce, but errors if happening not at the start of the simulation.
        """
        if self._built:
            raise ValueError("Virtual sites cannot be added during the simulation.")
        return super()._add_bond(members, params)

    @classmethod
    def uses_pbc(cls) -> bool:
        """Whether the force uses periodic boundary conditions. Returns false for all virtual sites."""
        return False

    def _remove_bond(self, bond_id: int) -> None:
        """Overridden by this class to throw errors, as removing virtual sites during simulations is not supported."""
        raise ValueError("Virtual sites cannot be removed during the simulation.")

    def _destroy(self) -> None:
        """Overridden by this class to throw errors, as removing virtual sites during simulations is not supported."""
        raise ValueError("Virtual sites can't be destoryed.")

    def _set_force_obj(self) -> mm.Force:
        """Unreachable function for virtual sites."""
        assert False  # unreachable

    @classmethod
    def _add_to_force(
        cls, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        """Unreachable function for virtual sites."""
        assert False  # unreachable
