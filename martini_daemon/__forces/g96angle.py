import openmm as mm
from ..__core import BondedForce, register_available_force
from ..__parser import register_angle_type


@register_angle_type(type_=2, args=["degree", "float"])
@register_available_force
class G96Angle(BondedForce):

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
        return "g96_angle"

    def _set_force_obj(self):
        self.force = mm.CustomAngleForce("0.5 * k * (cos(theta) - cos(theta0))^2")
        self.force.addPerAngleParameter("theta0")
        self.force.addPerAngleParameter("k")

    filters = {"angle"}
