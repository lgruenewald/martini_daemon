import openmm as mm

from ..__core import BondedForce, register_available_force
from ..__parser import register_angle_type


@register_angle_type(type_=1, args=["degree", "float"])
@register_available_force
class HarmonicAngle(BondedForce):
    @classmethod
    def _add_to_force(
        cls, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        assert isinstance(force, mm.HarmonicAngleForce)
        force.addAngle(*members, *params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @staticmethod
    def uses_pbc() -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "harmonic_angle"

    def _set_force_obj(self) -> mm.Force:
        return mm.HarmonicAngleForce()

    filters = {"angle"}
