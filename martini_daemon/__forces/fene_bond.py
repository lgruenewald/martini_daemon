import openmm as mm
from ..__core import BondedForce, register_available_force
from ..__parser import register_bond_type


@register_bond_type(type_=7, args=["float", "float"], is_excl=True)
@register_available_force
class FENEBond(BondedForce):
    # FENE (finitely extensible nonlinear elastic) bond
    def _add_to_force(self, members: list[int], params: list[float]) -> None:
        self.force.addBond(*members, params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @staticmethod
    def uses_pbc() -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "fene_bond"

    def _set_force_obj(self):
        self.force = mm.CustomBondForce("- 0.5 * k * b^2 * log(1 - r^2 / b^2)")
        self.force.addPerBondParameter("b")  # equilibrium length
        self.force.addPerBondParameter("k")  # force constant

    filters = {"bond"}
