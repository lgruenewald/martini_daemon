import openmm as mm

from ..__core import VirtualSite, register_available_force
from ..__parser import register_vsiten_type


@register_vsiten_type(3, ["float"])
@register_vsiten_type(1, [])
@register_available_force
class VSiteWeightedAverage(VirtualSite):
    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        n = len(members) - 1
        if len(params) == 0:
            return [1.0 / n] * n
        assert len(params) == n
        return [w / sum(params) for w in params]

    @classmethod
    def get_name(cls) -> str:
        return "weighted_average"

    def _make_vsite(self, vid, other, params) -> mm.VirtualSite:
        n = len(other)
        if n == 1:
            return mm.TwoParticleAverageSite(other[0], other[0], 1.0, 0.0)
        if n == 2:
            return mm.TwoParticleAverageSite(other[0], other[1], params[0], params[1])
        if n == 3:
            return mm.ThreeParticleAverageSite(
                other[0], other[1], other[2], params[0], params[1], params[2]
            )
        return mm.LocalCoordinatesSite(
            other,  # atoms
            params,  # origin weights
            [0.0] * n,  # x direction weight
            [0.0] * n,  # y direction weight
            [0.0, 0.0, 0.0],  # coordinates
        )

    filters = {"virtual_site", "vsite"}
