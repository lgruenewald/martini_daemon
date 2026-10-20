import openmm as mm

from ..__core import (
    BondedForce,
    System,
    register_available_force,
    sigma_epsilon_to_c6_c12,
)
from ..__parser import register_pair_type


@register_pair_type(type_=1, args=["float", "float"])
@register_available_force
class Pairs(BondedForce):
    def __init__(self, system: System) -> None:
        super().__init__(system)
        self.epsilon_r = system.additional_data["epsilon_r"]

    @classmethod
    def _add_to_force(
        cls, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        assert isinstance(force, mm.CustomBondForce)
        force.addBond(*members, params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        q1 = self.system.get_charge(members[0])
        q2 = self.system.get_charge(members[1])
        fudgeQQ = self.system.additional_data.get("fudgeQQ")
        q_prod = q1 * q2 * fudgeQQ

        if self.system.additional_data.get("combination_rule") == 1:
            c6, c12 = params
        else:
            c6, c12 = sigma_epsilon_to_c6_c12(*params)

        return [q_prod, c6, c12]

    @classmethod
    def uses_pbc(cls) -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "pair"

    def _set_force_obj(self) -> mm.Force:
        force = mm.CustomBondForce(
            "LJ + ES;"
            "LJ = (C12 / r^12 - C6 / r^6);"
            "ES = f*qprod/epsilon_r/r;"
            f"epsilon_r = {self.epsilon_r};"
            "f = 138.935458;"
        )
        force.addPerBondParameter("qprod")
        force.addPerBondParameter("C6")
        force.addPerBondParameter("C12")

        return force
