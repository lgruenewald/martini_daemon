import openmm as mm

from ..__core import BondedForce, register_available_force
from ..__parser import register_angle_type


@register_angle_type(type_=9, args=["float", "float"])
@register_available_force
class LinearAngle(BondedForce):
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
        return False

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "linear_angle"

    def _set_force_obj(self) -> mm.Force:
        # FIXME: this is not pbc friendly
        force = mm.CustomCompoundBondForce(
            3,
            "0.5*k*distj2; "
            "distj2=(xj-x2)^2+(yj-y2)^2+(zj-z2)^2; "
            "xj=a*x1+(1-a)*x3; "
            "yj=a*y1+(1-a)*y3; "
            "zj=a*z1+(1-a)*z3;",
        )
        force.addPerBondParameter("a")
        force.addPerBondParameter("k")

        return force

    filters = {"angle"}
