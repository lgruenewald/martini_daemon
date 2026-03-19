import openmm as mm
from ..__core import BondedForce, register_available_force
from ..__parser import register_angle_type

@register_angle_type(type_=9, args=["float", "float"])
@register_available_force
class LinearAngle(BondedForce):

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
        return "linear_angle"

    def _set_force_obj(self):
        self.force = mm.CustomCompoundBondForce(
            3,  # 3 particles per compund bond force
            "0.5*k*distj2; "
            "distj2=(xj-x2)^2+(yj-y2)^2+(zj-z2)^2; "
            "xj=a*x1+(1-a)*x3; "
            "yj=a*y1+(1-a)*y3; "
            "zj=a*z1+(1-a)*z3;"
        )
        self.force.addPerBondParameter("a")
        self.force.addPerBondParameter("k")

    filters = {"angle"}
