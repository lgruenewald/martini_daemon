import openmm as mm
from .vsite import VirtualSite


class VSiteThree(VirtualSite):

    def _make_vsite(self, vid, members, weights):
        vsite = mm.ThreeParticleAverageSite(
            members[0], members[1], members[2],
            weights[0], weights[1], weights[2]
        )
        self._sysstar._system.setVirtualSite(vid, vsite)

    def is_instance(self, filter):
        return filter in {"virtual_site", "vsite", "vsite3"}
