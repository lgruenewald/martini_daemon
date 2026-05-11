import openmm as mm

from ..__core import BondedForce, register_available_force
from ..__parser import register_bond_type


@register_bond_type(type_=10, args=["float", "float", "float", "float"], is_excl=True)
@register_available_force
class DistanceRestraint(BondedForce):
    @classmethod
    def _add_to_force(
        cls, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        assert isinstance(force, mm.CustomBondForce)
        force.addBond(*members, params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @classmethod
    def uses_pbc(cls) -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "distance_restraint"

    def _set_force_obj(self) -> mm.Force:
        force = mm.CustomBondForce(
            "V_low*step(low-r)+V_up1*step(r-up1)*step(up2-r)+V_up2*step(r-up2); "
            "V_low = 0.5 * k *(r-low)^2; "
            "V_up1 = 0.5 * k * (r-up1)^2; "
            "V_up2 = 0.5 * k * (up2 - up1) * (2*r-up2-up1);"
        )
        force.addPerBondParameter("low")
        force.addPerBondParameter("up1")
        force.addPerBondParameter("up2")
        force.addPerBondParameter("k")

        return force

    filters = {"bond"}
