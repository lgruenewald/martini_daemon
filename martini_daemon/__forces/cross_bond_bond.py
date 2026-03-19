import openmm as mm
from ..__core import BondedForce, register_available_force
from ..__parser import register_angle_type

@register_angle_type(type_=3, args=["float", "float", "float"])
@register_available_force
class CrossBondBond(BondedForce):
    def _add_to_force(self, members: list[int], params: list[float]) -> None:
        self.force.addBond(members, params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @staticmethod
    def uses_pbc() -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "cross_bond_bond"

    def _set_force_obj(self):
        self.force = mm.CustomCompoundBondForce(
            3,  # 3 particles per compund bond force
            "k*(distance(p1,p2)-r1)*(distance(p3,p2)-r2)"
        )
        self.force.addPerBondParameter("r1")
        self.force.addPerBondParameter("r2")
        self.force.addPerBondParameter("k")

    filters = {"angle"}
