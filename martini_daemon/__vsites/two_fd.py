import openmm as mm
from ..__core import VirtualSite, register_available_force
from ..__parser import register_vsite2_type


@register_vsite2_type(2, ["float"])
@register_available_force
class VSite2fd(VirtualSite):
    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @classmethod
    def get_name(cls) -> str:
        return "2fd"

    def _make_vsite(self, vid, members, params):
        i, j = members
        return mm.LocalCoordinatesSite(
            [i, j, j],  # atoms
            [1.0, 0.0, 0.0],  # origin weights
            [-0.5, 0.5, 0.0],  # x direction weight
            [0.0, 0.0, 0.0],  # y direction weight
            [params[0], 0.0, 0.0],  # coordinatesI
        )

    filters = {"virtual_site", "vsite"}
