from ..__core import BondedForce, register_available_force
from ..__parser import register_bond_type

@register_bond_type(type_=5, args=[], is_excl=True)
@register_available_force
class Connection(BondedForce):
    def build(self, must=False) -> None:
        pass

    def _set_force_obj(self) -> None:
        pass

    def _add_to_force(self, members: list[int], params: list[float]) -> None:
        pass

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @staticmethod
    def uses_pbc() -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "connection"

    def build(self):
        pass

    def _destroy(self):
        pass

    filters = {"bond"}
