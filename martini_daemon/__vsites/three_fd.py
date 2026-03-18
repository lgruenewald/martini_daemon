import openmm as mm
from ..__core import VirtualSite

class VSite3fd(VirtualSite):
    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @classmethod
    def get_name(cls) -> str:
        return "3fd"

    def _make_vsite(self, vid, members, params):
        a, d = params
        return mm.LocalCoordinatesSite(
            members,  # atoms
            [1.0, 0.0, 0.0],  # origin weights
            [-1.0, 1.0 - a, a],  # x direction weight
            [0.0, 0.0, 0.0],  # y direction weight
            [d, 0.0, 0.0]  # coordinates
        )

    filters = {"virtual_site", "vsite"}
