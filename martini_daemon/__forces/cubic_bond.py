import openmm as mm
from ..__core import BondedForce, register_available_force
from ..__parser import register_bond_type


@register_bond_type(type_=4, args=["float", "float", "float"], is_excl=True)
@register_available_force
class CubicBond(BondedForce):
    def _add_to_force(self, members: list[int], params: list[float]) -> None:
        self.force.addBond(*members, params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @staticmethod
    def uses_pbc() -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    def _set_force_obj(self) -> None:
        self.force = mm.CustomBondForce("kb * (r - b)^2 + kb * kcub * (r - b)^3")
        self.force.addPerBondParameter("b")  # equilibrium length
        self.force.addPerBondParameter("kb")  # force constant
        self.force.addPerBondParameter("kcub")  # cubic force constant

    @classmethod
    def get_name(cls) -> str:
        return "cubic_bond"

    filters = {"bond"}
