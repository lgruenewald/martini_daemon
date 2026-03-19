import openmm as mm
from ..__core import VirtualSite, register_available_force
from ..__parser import register_vsite3_type

@register_vsite3_type(1, ["float", "float"])
@register_available_force
class VSiteThree(VirtualSite):

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        a, b = params
        return [1. - a - b, a, b]

    @classmethod
    def get_name(cls) -> str:
        return "vsite3"

    def _make_vsite(self, vid, members, weights):
        return mm.ThreeParticleAverageSite(
            members[0], members[1], members[2],
            weights[0], weights[1], weights[2]
        )

    filters = {"virtual_site", "vsite"}
