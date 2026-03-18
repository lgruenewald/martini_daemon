import openmm as mm
from ..__core import VirtualSite

class VSiteCenterOfMass(VirtualSite):

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @classmethod
    def get_name(cls) -> str:
        return "center_of_mass"

    def _make_vsite(self, vid, others, params) -> mm.VirtualSite:
        n = len(others)
        masses = [
            self.system.get_mass(i)
            for i in others
        ]
        weights = [m/sum(masses) for m in masses]
        return mm.LocalCoordinatesSite(
            others,  # atoms
            weights,  # origin weights
            [0.0] * n,  # x direction weight
            [0.0] * n,  # y direction weight
            [0.0, 0.0, 0.0]  # coordinates
        )

    filters = {"virtual_site", "vsite", "com"}
