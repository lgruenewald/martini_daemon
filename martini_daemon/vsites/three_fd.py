import openmm as mm
from .vsite import VirtualSite


class VSite3fd(VirtualSite):
    def _make_vsite(self, vid, members, params):
        i, j, k = members
        a, d = params
        vsite = mm.LocalCoordinatesSite(
            [i, j, k],  # atoms
            [1.0, 0.0, 0.0],  # origin weights
            [-1.0, 1.0 - a, a],  # x direction weight
            [0.0, 0.0, 0.0],  # y direction weight
            [d, 0.0, 0.0]  # coordinates
        )
        self._sysstar._system.setVirtualSite(vid, vsite)

    _filters = {"virtual_site", "vsite", "3fd"}
