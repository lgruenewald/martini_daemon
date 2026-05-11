import openmm as mm

from ..__core import VirtualSite, register_available_force
from ..__parser import register_vsite4_type


@register_vsite4_type(2, ["float", "float", "float"])
@register_available_force
class VSite4fdn(VirtualSite):
    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @classmethod
    def get_name(cls) -> str:
        return "4fdn"

    def _make_vsite(
        self, vid: int, other: list[int], params: list[float]
    ) -> mm.VirtualSite:
        a, b, c = params
        return mm.LocalCoordinatesSite(
            other,  # atoms
            [1.0, 0.0, 0.0, 0.0],  # origin weights
            [1.0 - a, -1.0, a, 0.0],  # x direction weight
            [1.0 - b, -1.0, 0.0, b],  # y direction weight
            [0.0, 0.0, c],  # coordinates
        )

    filters = {"virtual_site", "vsite"}
