import openmm as mm

from ..__core import BondedForce, register_available_force
from ..__parser import register_dihedral_type


@register_dihedral_type(type_=3, args=["float" for _ in range(6)])
@register_available_force
class RBTorsion(BondedForce):
    @classmethod
    def _add_to_force(
        cls, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        assert isinstance(force, mm.RBTorsionForce)
        force.addTorsion(*members, *params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @staticmethod
    def uses_pbc() -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "rb_torsion"

    def _set_force_obj(self) -> mm.Force:
        return mm.RBTorsionForce()

    filters = {"dihedral"}


@register_dihedral_type(type_=5, args=["float" for _ in range(4)])
@register_available_force
class FourierDihedral(BondedForce):
    @classmethod
    def _add_to_force(
        cls, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        assert isinstance(force, mm.RBTorsionForce)
        force.addTorsion(*members, *params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return [
            params[1] + 0.5 * (params[0] + params[2]),
            0.5 * (-params[0] + 3 * params[2]),
            -params[1] + 4 * params[3],
            -2 * params[2],
            -4 * params[3],
            0.0,
        ]

    @staticmethod
    def uses_pbc() -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "fourier_dihedral"

    def _set_force_obj(self) -> mm.Force:
        return mm.RBTorsionForce()

    filters = {"dihedral"}
