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

    def _make_vsite(
        self, vid: int, other: list[int], params: list[float]
    ) -> mm.VirtualSite:
        return mm.TwoParticleAverageSite(other[0], other[0], 1.0, 0.0)

    filters = {"virtual_site", "vsite"}
