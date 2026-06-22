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
        if type_int == 1:
            return "posres"
        if type_int == 2:
            return "flat_bottomed_restraint"
        raise ParseException(f"Unsupported position restraint type {type_int}.")

    @classmethod
    def get_type_args(cls, type_int: int) -> list[str]:
        if type_int == 1:
            return ["float", "float", "float"]
        if type_int == 2:
            return ["int", "float", "float"]
        raise ParseException(f"Unsupported position restraint type {type_int}.")

    @classmethod
    def get_name(cls) -> str:
        return "position_restraints"


@register_available_force
class FlatBottomedRestraint(BondedForce):
    @classmethod
    def _add_to_force(
        cls, force: mm.Force, members: list[int], params: list[float]
    ) -> None:
        assert isinstance(force, mm.CustomExternalForce)
        assert len(members) == 1
        i = members[0]
        force.addParticle(i, params)        

    @classmethod
    def uses_pbc(cls) -> bool:
        return False

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "flat_bottomed_restraint"

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        assert len(members) == 1
        i = members[0]
        x0, y0, z0 = self.system.additional_data["respos"][i]
        g, r, k = params
        if g == 2:
            g = 8
        if g not in {1, 3, 4, 5, 6, 7, 8}:
            raise ValueError(
                "The g parameter for flat bottomed position restraints must be in {1, 3, 4, 5, 6, 7, 8}."
            )
        rsign = 1
        if r < 0:
            r *= -1
            rsign = -1
        return [g, r, k, rsign, x0, y0, z0]

    def _set_force_obj(self) -> mm.Force:
        force = mm.CustomExternalForce(
            "0.5*k*(d-r)^2*step(rsign*(d-r));"
            "d=delta(g-1)*sphere+delta(g-3)*layerx+delta(g-4)*layery+delta(g-5)*layerz+delta(g-6)*cylx+delta(g-7)*cyly+delta(g-8)*cylz;"
            "sphere=periodicdistance(x,y,z,x0,y0,z0);"
            "layerx=periodicdistance(x,0,0,x0,0,0);"
            "layery=periodicdistance(0,y,0,0,y0,0);"
            "layerz=periodicdistance(0,0,z,0,0,z0);"
            "cylx=sqrt(periodicdistance(0,y,0,0,y0,0)^2+periodicdistance(0,0,z,0,0,z0)^2);"
            "cyly=sqrt(periodicdistance(x,0,0,x0,0,0)^2+periodicdistance(0,0,z,0,0,z0)^2);"
            "cylz=sqrt(periodicdistance(0,y,0,0,y0,0)^2+periodicdistance(x,0,0,x0,0,0)^2);"
        )
        force.addPerParticleParameter("g")
        force.addPerParticleParameter("r")
        force.addPerParticleParameter("k")
        force.addPerParticleParameter("rsign")
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")

        return force


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

    @classmethod
    def uses_pbc(cls) -> bool:
        # FIXME: WHY DOES THIS WORK
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
