import openmm as mm

from ..__core import BondedForce, register_available_force
from ..__parser import register_bond_type


@register_bond_type(type_=1, args=["float", "float"], is_excl=True)
@register_bond_type(type_=6, args=["float", "float"], is_excl=False)
@register_available_force
class HarmonicBond(BondedForce):
    def _add_to_force(
        self, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        assert isinstance(force, mm.HarmonicBondForce)
        force.addBond(*members, *params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @staticmethod
    def uses_pbc() -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    def _set_force_obj(self) -> None:
        self.force = mm.HarmonicBondForce()

    filters = {"bond"}

    @classmethod
    def get_name(cls) -> str:
        return "harmonic_bond"
