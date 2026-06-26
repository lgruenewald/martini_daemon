import openmm as mm

from ..__core import BondedForce, register_available_force
from ..__parser import register_angle_type


@register_angle_type(type_=3, args=["float", "float", "float"])
@register_available_force
class CrossBondBond(BondedForce):
    @classmethod
    def _add_to_force(
        cls, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        assert isinstance(force, mm.CustomCompoundBondForce)
        force.addBond(members, params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @classmethod
    def uses_pbc(cls) -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "cross_bond_bond"

    def _set_force_obj(self) -> mm.Force:
        force = mm.CustomCompoundBondForce(
            3,
            "k*(distance(p1,p2)-r1)*(distance(p3,p2)-r2)",
        )
        force.addPerBondParameter("r1")
        force.addPerBondParameter("r2")
        force.addPerBondParameter("k")
        return force

    filters = {"angle"}
