import openmm as mm

from ..__core import VirtualSite, register_available_force
from ..__parser import register_vsite3_type


@register_vsite3_type(4, ["float", "float", "float"])
@register_available_force
class VSite3out(VirtualSite):
    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @classmethod
    def get_name(cls) -> str:
        return "3out"

    def _make_vsite(
        self, vid: int, other: list[int], params: list[float]
    ) -> mm.VirtualSite:
        i, j, k = other
        a, b, c = params
        return mm.OutOfPlaneSite(i, j, k, a, b, c)

    filters = {"virtual_site", "vsite"}
