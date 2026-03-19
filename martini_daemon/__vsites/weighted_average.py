import openmm as mm

from ..__parser import register_vsiten_type
from ..__core import VirtualSite, register_available_force

@register_vsiten_type(3, ["float"])
@register_vsiten_type(1, [])
@register_available_force
class VSiteWeightedAverage(VirtualSite):

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        n = len(members) - 1
        if len(params) == 0:
            return [1./n for _ in range(n)]
        else:
            assert len(params) == n
            return [w/sum(params) for w in params]

    @classmethod
    def get_name(cls) -> str:
        return "weighted_average"

    def _make_vsite(self, vid, members, weights):
        n = len(members)
        if n == 1:
            return mm.TwoParticleAverageSite(
                members[0], members[0], 1.0, 0.0
            )
        elif n == 2:
            return mm.TwoParticleAverageSite(
                members[0], members[1], weights[0], weights[1]
            )
        elif n == 3:
            return mm.ThreeParticleAverageSite(
                members[0], members[1], members[2],
                weights[0], weights[1], weights[2]
            )
        else:
            return mm.LocalCoordinatesSite(
                members,  # atoms
                weights,  # origin weights
                [0.0] * n,  # x direction weight
                [0.0] * n,  # y direction weight
                [0.0, 0.0, 0.0]  # coordinates
            )

    filters = {"virtual_site", "vsite"}
