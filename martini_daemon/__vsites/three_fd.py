import openmm as mm

from ..__core import VirtualSite, register_available_force
from ..__parser import register_vsite3_type


@register_vsite3_type(2, ["float", "float"])
@register_available_force
class VSite3fd(VirtualSite):
    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @classmethod
    def get_name(cls) -> str:
        return "3fd"

    def _make_vsite(self, vid, other, params) -> mm.VirtualSite:
        a, d = params
        return mm.LocalCoordinatesSite(
            other,  # atoms
            [1.0, 0.0, 0.0],  # origin weights
            [-1.0, 1.0 - a, a],  # x direction weight
            [0.0, 0.0, 0.0],  # y direction weight
            [d, 0.0, 0.0],  # coordinates
        )

    filters = {"virtual_site", "vsite"}
