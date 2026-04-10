import openmm as mm

from ..__core import BondedForce, register_available_force
from ..__parser import register_dihedral_type


@register_dihedral_type(type_=11, args=["float"] * 6)
@register_available_force
class CombinedBendingTorsion(BondedForce):
    def _add_to_force(
        self, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        assert isinstance(force, mm.CustomCompoundBondForce)
        force.addBond(members, params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @staticmethod
    def uses_pbc() -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "combined_bending_torsion"

    def _set_force_obj(self):
        self.force = mm.CustomCompoundBondForce(
            4,
            "k*sintheta0^3*sintheta1^3*(a0 + a1*cosphi + a2*cosphi^2 + a3*cosphi^3 + a4*cosphi^4); "
            "sintheta0 = sin(angle(p1, p2, p3));"
            "sintheta1 = sin(angle(p2, p3, p4));"
            "cosphi = cos(dihedral(p1, p2, p3, p4));",
        )
        self.force.addPerBondParameter("k")
        self.force.addPerBondParameter("a0")
        self.force.addPerBondParameter("a1")
        self.force.addPerBondParameter("a2")
        self.force.addPerBondParameter("a3")
        self.force.addPerBondParameter("a4")

    filters = {"dihedral"}
