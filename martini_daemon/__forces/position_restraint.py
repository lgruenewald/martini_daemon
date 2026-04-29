import openmm as mm

from ..__core import BondedForce, register_available_force
from ..__parser import InteractionDirective, ParseException, register_directive


@register_directive
class PositionRestraintDirective(InteractionDirective):
    @classmethod
    def get_number_members(cls) -> int:
        return 1

    @classmethod
    def get_number_params(cls, type_: int) -> tuple[int, int]:
        return 3, 3

    @classmethod
    def get_type(cls, type_int: int) -> str | None:
        if type_int != 1:
            raise ParseException(f"Unsupported position restraint type {type_int}.")
        return "posres"

    @classmethod
    def get_type_args(cls, type_int: int) -> list[str]:
        return ["float", "float", "float"]

    @classmethod
    def get_name(cls) -> str:
        return "position_restraints"


@register_available_force
class PositionRestraint(BondedForce):
    @classmethod
    def _add_to_force(
        cls, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        assert isinstance(force, mm.CustomExternalForce)
        assert len(members) == 1
        i = members[0]
        force.addParticle(i, params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        assert len(members) == 1
        i = members[0]
        kx, ky, kz = params
        x0, y0, z0 = self.system.additional_data["respos"][i]
        return [kx, ky, kz, x0, y0, z0]

    @staticmethod
    def uses_pbc() -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "posres"

    def _set_force_obj(self) -> mm.Force:
        # FIXME: this only works for 90 degree angles and when it doesn't grow/shrink
        force = mm.CustomExternalForce(
            "0.5*(fx+fy+fz);"
            "fx=kx*periodicdistance(x,0,0,x0,0,0)^2;"
            "fy=ky*periodicdistance(0,y,0,0,y0,0)^2;"
            "fz=kz*periodicdistance(0,0,z,0,0,z0)^2;"
        )
        force.addPerParticleParameter("kx")
        force.addPerParticleParameter("ky")
        force.addPerParticleParameter("kz")
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")

        return force
