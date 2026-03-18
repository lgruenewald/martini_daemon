import openmm as mm
from ..__core import VirtualSite


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
