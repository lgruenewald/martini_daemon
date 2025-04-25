import openmm as mm
from .vsite import VirtualSite


class VSiteOne(VirtualSite):

    def _make_vsite(self, vid, members, _):
        vsite = mm.TwoParticleAverageSite(
            members[0], members[0], 1.0, 0.0
        )
        self._sysstar._system.setVirtualSite(vid, vsite)

    def is_instance(self, filter):
        return filter in {"virtual_site", "vsite", "vsite1"}
