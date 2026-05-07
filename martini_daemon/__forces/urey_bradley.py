import openmm as mm

from ..__core import BondedForce, register_available_force
from ..__parser import register_angle_type


@register_angle_type(type_=5, args=["degree", "float", "float", "float"])
@register_available_force
class UreyBradley(BondedForce):
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
        return "urey_bradley"

    def _set_force_obj(self) -> mm.Force:
        force = mm.CustomCompoundBondForce(
            3,  # 3 particles per compund bond force
            "0.5*k*((angle(p1,p2,p3)-theta0)^2)+"  # angle part
            "0.5*kUB*((distance(p1,p3)-r13)^2)",  # distance part
        )
        force.addPerBondParameter("theta0")
        force.addPerBondParameter("k")
        force.addPerBondParameter("r13")
        force.addPerBondParameter("kUB")

        return force

    filters = {"angle"}
