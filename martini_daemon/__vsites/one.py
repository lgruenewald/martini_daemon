import openmm as mm
from ..__core import VirtualSite, register_available_force
from ..__parser import register_vsite1_type

@register_vsite1_type(1, [])
@register_available_force
class VSiteOne(VirtualSite):

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @classmethod
    def get_name(cls) -> str:
        return "vsite1"

    def _make_vsite(self, vid, members, _):
        return mm.TwoParticleAverageSite(
            members[0], members[0], 1.0, 0.0
        )

    filters = {"virtual_site", "vsite"}
