import openmm as mm

from ..__core import BondedForce, register_available_force
from ..__parser import register_bond_type


@register_bond_type(type_=5, args=[], is_excl=True)
@register_available_force
class Connection(BondedForce):
    def build(self, must: bool = False) -> None:
        pass

    def _set_force_obj(self) -> mm.Force:
        assert False

    @classmethod
    def _add_to_force(
        cls, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        pass

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @classmethod
    def uses_pbc(cls) -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "connection"

    def _destroy(self) -> None:
        pass

    filters = {"bond"}
