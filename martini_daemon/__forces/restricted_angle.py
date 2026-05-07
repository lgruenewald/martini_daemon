import openmm as mm

from ..__core import BondedForce, register_available_force
from ..__parser import register_angle_type


@register_angle_type(type_=10, args=["degree", "float"])
@register_available_force
class RestrictedAngle(BondedForce):
    @classmethod
    def _add_to_force(
        cls, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        assert isinstance(force, mm.CustomAngleForce)
        force.addAngle(*members, params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @classmethod
    def uses_pbc(cls) -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "restricted_angle"

    def _set_force_obj(self) -> mm.Force:
        force = mm.CustomAngleForce(
            "0.5*k*(cos(theta)-cos(theta0))^2/sin(theta)^2"
        )
        force.addPerAngleParameter("theta0")
        force.addPerAngleParameter("k")

        return force

    filters = {"angle"}
