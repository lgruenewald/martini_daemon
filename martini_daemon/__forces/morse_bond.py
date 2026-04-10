import openmm as mm

from ..__core import BondedForce, register_available_force
from ..__parser import register_bond_type


@register_bond_type(type_=3, args=["float", "float", "float"], is_excl=True)
@register_available_force
class MorseBond(BondedForce):
    def _add_to_force(
        self, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        assert isinstance(force, mm.CustomBondForce)
        force.addBond(*members, params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @staticmethod
    def uses_pbc() -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    def _set_force_obj(self) -> None:
        self.force = mm.CustomBondForce("D * (1 - exp(-beta * (r - b)))^2")
        self.force.addPerBondParameter("b")  # equilibrium length
        self.force.addPerBondParameter("D")  # force constant
        self.force.addPerBondParameter("beta")  # cubic force constant

    @classmethod
    def get_name(cls) -> str:
        return "morse_bond"

    filters = {"bond"}
