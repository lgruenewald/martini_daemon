import openmm as mm
from .vsite import VirtualSite


class VSite4fdn(VirtualSite):
    def _make_vsite(self, vid, members, params):
        i, j, k, L = members
        a, b, c = params
        vsite = mm.LocalCoordinatesSite(
            [i, j, k, L],  # atoms
            [1.0, 0.0, 0.0, 0.0],  # origin weights
            [1.0 - a, -1.0, a, 0.0],  # x direction weight
            [1.0 - b, -1.0, 0.0, b],  # y direction weight
            [0.0, 0.0, c]  # coordinates
        )
        self._sysstar._system.setVirtualSite(vid, vsite)

    def is_instance(self, filter):
        return filter in {"virtual_site", "vsite", "4fdn"}
