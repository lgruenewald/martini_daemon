import openmm as mm
from ..__core import BondedForce, register_available_force
from ..__parser import register_angle_type


@register_angle_type(type_=5, args=["degree", "float", "float", "float"])
@register_available_force
class UreyBradley(BondedForce):
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
        return "urey_bradley"

    def _set_force_obj(self):
        self.force = mm.CustomCompoundBondForce(
            3,  # 3 particles per compund bond force
            "0.5*k*((angle(p1,p2,p3)-theta0)^2)+"  # angle part
            "0.5*kUB*((distance(p1,p3)-r13)^2)",  # distance part
        )
        self.force.addPerBondParameter("theta0")
        self.force.addPerBondParameter("k")
        self.force.addPerBondParameter("r13")
        self.force.addPerBondParameter("kUB")

    filters = {"angle"}
