from math import cos, sin

import openmm as mm

from ..__core import VirtualSite, register_available_force
from ..__parser import register_vsite3_type


@register_vsite3_type(3, ["degree", "float"])
@register_available_force
class VSite3fad(VirtualSite):
    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @classmethod
    def get_name(cls) -> str:
        return "3fad"

    def _make_vsite(self, vid, other, params) -> mm.VirtualSite:
        theta, d = params
        return mm.LocalCoordinatesSite(
            other,  # atoms
            [1.0, 0.0, 0.0],  # origin weights
            [-0.5, 0.5, 0.0],  # x direction weight
            [0.0, -0.5, 0.5],  # y direction weight
            [d * cos(theta), d * sin(theta), 0.0],  # coordinates
        )

    filters = {"virtual_site", "vsite"}
