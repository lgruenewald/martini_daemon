from .bonded_force import BondedForce
from abc import ABCMeta, abstractmethod
import openmm as mm


class VirtualSite(BondedForce, metaclass=ABCMeta):

    def __init__(self, system):
        super().__init__(system)
        self._built = False

    @abstractmethod
    def _make_vsite(self, vid, other, params) -> mm.VirtualSite:
        raise NotImplementedError

    def build(self, must=False) -> None:
        if self._built:
            return
        self._built = True
        for index, (members, params) in self.iterate_bonds():
            vid, *other = members
            self.system._add_vsite(vid, self._make_vsite(vid, other, params))

    def delta_degrees_of_freedom(self) -> int:
        return self.num_bonds() * 3

    def _add_bond(self, members: list[int], params: list[float]) -> int:
        if self._built:
            raise ValueError("Virtual sites cannot be added during the simulation.")
        return super()._add_bond(members, params)

    @staticmethod
    def uses_pbc() -> bool:
        return False

    def _remove_bond(self, bond_id: int) -> None:
        raise ValueError("Virtual sites cannot be removed during the simulation.")

    def _destroy(self) -> None:
        raise ValueError("Virtual sites can't be destoryed.")

    def _set_force_obj(self) -> None:
        assert False  # unreachable

    def _add_to_force(self, members: list[int], params: list[float]) -> None:
        assert False  # unreachable
