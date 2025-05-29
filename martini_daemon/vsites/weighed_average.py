import openmm as mm
from .vsite import VirtualSite


class VSiteWeighedAverage(VirtualSite):

    def _make_vsite(self, vid, members, weights):
        n = len(members)
        if n == 1:
            vsite = mm.TwoParticleAverageSite(
                members[0], members[0], 1.0, 0.0
            )
        elif n == 2:
            vsite = mm.TwoParticleAverageSite(
                members[0], members[1], weights[0], weights[1]
            )
        elif n == 3:
            vsite = mm.ThreeParticleAverageSite(
                members[0], members[1], members[2],
                weights[0], weights[1], weights[2]
            )
        else:
            vsite = mm.LocalCoordinatesSite(
                members,  # atoms
                weights,  # origin weights
                [0.0] * n,  # x direction weight
                [0.0] * n,  # y direction weight
                [0.0, 0.0, 0.0]  # coordinates
            )
        self._sysstar._system.setVirtualSite(vid, vsite)

    def is_instance(self, filter):
        return filter in {"virtual_site", "vsite", "weighed_average"}
