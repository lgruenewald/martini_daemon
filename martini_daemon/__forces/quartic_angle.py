import openmm as mm
from ..__core import BondedForce, register_available_force
from ..__parser import register_angle_type

@register_angle_type(type_=6, args=["degree", "float", "float", "float", "float", "float"])
@register_available_force
class QuarticAngle(BondedForce):
    def _add_to_force(self, members: list[int], params: list[float]) -> None:
        self.force.addAngle(*members, params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @staticmethod
    def uses_pbc() -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "quartic_angle"

    def _set_force_obj(self):
        self.force = mm.CustomAngleForce(
            "c0+"
            "c1*(theta-theta0)+"
            "c2*(theta-theta0)^2+"
            "c3*(theta-theta0)^3+"
            "c4*(theta-theta0)^4"
        )
        self.force.addPerAngleParameter("theta0")
        self.force.addPerAngleParameter("c0")
        self.force.addPerAngleParameter("c1")
        self.force.addPerAngleParameter("c2")
        self.force.addPerAngleParameter("c3")
        self.force.addPerAngleParameter("c4")

    filters = {"angle"}
