from __future__ import annotations
from abc import ABCMeta, abstractmethod
from typing import Iterable
from .force import Force

class BondedForce(Force, metaclass=ABCMeta):

    def __init__(self, system):
        super().__init__(system)
        self.entries: dict[int, tuple[list[int], list[float]]] = {}

    def prepare_force_obj(self) -> None:
        super().prepare_force_obj()
        if self.pbc:
            self.force.setUsesPeriodicBoundaryConditions(True)
        # TODO add all from entries

    @abstractmethod
    def parse(self, members: list[int], params: list[float]) -> list[float]:
        """
        Given a list of members and parameters (as unwrapped = minimal pre-parsing), pa
        """
        raise NotImplementedError

    @abstractmethod
    def passes_filter(self, filter_: str) -> bool:
        return filter_ in self.filters

    # list of categories in which this bonded force should be included, including broad categories e.g. "bond"
    filters: set[str] = set()
    # abstract fields:
    members: int
    pbc: bool

    def add_bond(self, members: list[int], params: list[float]) -> int:
        raise NotImplementedError

    def remove_bond(self, bond_id: int) -> None:
        raise NotImplementedError

    def iterate_bonds(self) -> Iterable[tuple[int, tuple[list[int], list[float]]]]:
        return self.entries.items()