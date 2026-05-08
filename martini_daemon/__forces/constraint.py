from typing import NoReturn

import openmm as mm

from ..__core import BondedForce, System, register_available_force
from ..__parser import register_constraint_type


@register_constraint_type(type_=1, args=["float"], is_excl=True)
@register_constraint_type(type_=2, args=["float"], is_excl=False)
@register_available_force
class Constraint(BondedForce):
    def __init__(self, system: System) -> None:
        super().__init__(system)
        self.__built = False

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    def build(self, must: bool = False) -> None:
        if not self.__built or must:
            self.__built = True
            for index, (members, params) in self.iterate_bonds():
                self.system._add_constraint(*members, *params)

    @classmethod
    def uses_pbc(cls) -> bool:
        return False

    def delta_degrees_of_freedom(self) -> int:
        # works for *most* constraint constructs I think
        return self.num_bonds()

    @classmethod
    def get_name(cls) -> str:
        return "constraint"

    filters = {"bond"}

    def _set_force_obj(self) -> mm.Force:
        assert False

    def remove(self, i: int) -> NoReturn:
        raise Exception("Can't remove constraints")

    def should_build(self) -> bool:
        assert False

    @classmethod
    def _add_to_force(
        cls, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        assert False

    def _destroy(self) -> NoReturn:
        raise Exception("Can't destroy constraints.")
