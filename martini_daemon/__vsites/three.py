import openmm as mm

from ..__core import VirtualSite, register_available_force
from ..__parser import register_vsite3_type


@register_vsite3_type(1, ["float", "float"])
@register_available_force
class VSiteThree(VirtualSite):
    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        a, b = params
        return [1.0 - a - b, a, b]

    @classmethod
    def get_name(cls) -> str:
        return "vsite3"

    def _make_vsite(self, vid, other, params) -> mm.VirtualSite:
        return mm.ThreeParticleAverageSite(
            other[0], other[1], other[2], params[0], params[1], params[2]
        )

    filters = {"virtual_site", "vsite"}
