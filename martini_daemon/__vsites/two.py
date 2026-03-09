import openmm as mm
from .vsite import VirtualSite


class VSiteTwo(VirtualSite):

    def _make_vsite(self, vid, members, weights):
        vsite = mm.TwoParticleAverageSite(
            members[0], members[1], weights[0], weights[1]
        )
        self._sysstar._system.setVirtualSite(vid, vsite)

    _filters = {"virtual_site", "vsite", "vsite2"}
