from __future__ import annotations
from abc import ABCMeta, abstractmethod
from typing import Iterable
from .force import Force

class BondedForce(Force, metaclass=ABCMeta):

    def __init__(self, system):
        super().__init__(system)
        self.__entries: dict[int, tuple[list[int], list[float]]] = {}
        self.__next_entry_id = 0

    def _prepare_force_obj(self) -> None:
        super()._prepare_force_obj()
        if self.uses_pbc():
            self.force.setUsesPeriodicBoundaryConditions(True)
        for members, params in self.__entries.values():
            self._add_to_force(members, params)

    @abstractmethod
    def _add_to_force(self, members: list[int], params: list[float]) -> None:
        """
        Protected because only this class should call this.
        """
        raise NotImplementedError

    def should_build(self) -> bool:
        """
        If True, build() will be called on a global rebuild.

        By default this happens if there are any entries and there is no force.

        To trigger a rebuild, you generally therefore want to call _destroy().
        """
        return len(self.__entries) > 0 and self.force is None

    @abstractmethod
    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        """
        Protected because only this class should call this.

        Given a list of members and parameters (as unwrapped = minimal pre-parsing), parse the params to how they
        should be stored in entries and passed to _add_to_force.
        """
        raise NotImplementedError

    def passes_filter(self, filter_: str) -> bool:
        return filter_ == self.get_name() or filter_ in self.filters

    # List of categories in which this bonded force should be included, including broad categories e.g. "bond"
    filters: set[str] = set()

    @staticmethod
    @abstractmethod
    def uses_pbc() -> bool:
        """
        Should return False if and only if this interaction is unable to handle being its constituent atoms be
        in different instances of the periodic box.

        Note: will call setUsesPeriodicBoundaryConditions on self.force. Override _prepare_force_obj for really
        custom force objects or no force objects that obey PBC but have no such method.
        """
        raise NotImplementedError

    def _add_bond(self, members: list[int], params: list[float]) -> int:
        """
        API function to parse parameters and add it to the force.

        Protected because only System should call this.
        """
        params = self._parse(members, params)
        self.__entries[self.__next_entry_id] = (members, params)
        if self.force is not None:
            self._add_to_force(members, params)
            self.system.flag_reinitialize()
        self.__next_entry_id += 1
        return self.__next_entry_id - 1

    def _remove_bond(self, bond_id: int) -> None:
        """
        Protected because only System should call this.
        """
        del self.__entries[bond_id]
        if self.force is not None:
            self._destroy()

    def iterate_bonds(self) -> Iterable[tuple[int, tuple[list[int], list[float]]]]:
        """
        Returns an Iterable over bond_id's and a tuple of members and params.

        WARNING! Do not edit the bonds given by this function!
        """
        return self.__entries.items()

    def num_bonds(self) -> int:
        return len(self.__entries)

    def get_members(self, bond_id: int) -> list[int]:
        members, _ = self.__entries[bond_id]
        return members

    @classmethod
    def is_coupling(cls):
        return False
