import openmm as mm
from ..__core import VirtualSite, register_available_force
from ..__parser import register_vsite2_type


@register_vsite2_type(1, ["float"])
@register_available_force
class VSiteTwo(VirtualSite):

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return [1 - params[0], params[0]]

    @classmethod
    def get_name(cls) -> str:
        return "vsite2"

    def _make_vsite(self, vid, members, weights):
        return mm.TwoParticleAverageSite(members[0], members[1], weights[0], weights[1])

    _filters = {"virtual_site", "vsite"}
