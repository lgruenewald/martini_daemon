import openmm as mm
from .vsite import VirtualSite


class VSite3out(VirtualSite):
    def _make_vsite(self, vid, members, params):
        i, j, k = members
        a, b, c = params
        vsite = mm.OutOfPlaneSite(
            i, j, k, a, b, c
        )
        self._sysstar._system.setVirtualSite(vid, vsite)

    _filters = {"virtual_site", "vsite", "3out"}
