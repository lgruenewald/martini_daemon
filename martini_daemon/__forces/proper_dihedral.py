import openmm as mm

from ..__core import BondedForce, register_available_force
from ..__parser import register_dihedral_type


@register_dihedral_type(type_=1, args=["degree", "float", "float"])
@register_dihedral_type(type_=4, args=["degree", "float", "float"])
@register_dihedral_type(type_=9, args=["degree", "float", "float"])
@register_available_force
class ProperDihedral(BondedForce):
    @classmethod
    def _add_to_force(
        cls, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        assert isinstance(force, mm.PeriodicTorsionForce)
        force.addTorsion(*members, *params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        theta, force, mult = params
        return [mult, theta, force]

    @classmethod
    def uses_pbc(cls) -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "proper_dihedral"

    def _set_force_obj(self) -> mm.Force:
        return mm.PeriodicTorsionForce()

    filters = {"dihedral"}
