import openmm as mm
import math
from ..__core import BondedForce, register_available_force
from ..__parser import register_dihedral_type

@register_dihedral_type(type_=2, args=["degree", "float"])
@register_available_force
class ImproperDihedral(BondedForce):

    def _add_to_force(self, members: list[int], params: list[float]) -> None:
        self.force.addTorsion(*members, params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @staticmethod
    def uses_pbc() -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "improper_dihedral"

    def _set_force_obj(self):
        self.force = mm.CustomTorsionForce(
            "0.5*k*(thetap-theta0)^2;"
            "thetap = step(-plus)*2*pi+theta+step(minus)*(-2*pi);"
            "plus=theta+pi-theta0;"
            "minus=theta-pi-theta0;"
            f"pi = {math.pi:.14f}"
        )
        self.force.addPerTorsionParameter("theta0")
        self.force.addPerTorsionParameter("k")


    _filters = {"dihedral"}
