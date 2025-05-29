import openmm as mm
from .vsite import VirtualSite


class VSite2fd(VirtualSite):
    def _make_vsite(self, vid, members, params):
        i, j = members
        (d,) = params
        vsite = mm.LocalCoordinatesSite(
            [i, j, j],  # atoms
            [1.0, 0.0, 0.0],  # origin weights
            [-0.5, 0.5, 0.0],  # x direction weight
            [0.0, 0.0, 0.0],  # y direction weight
            [d, 0.0, 0.0]  # coordinatesI
        )
        self._sysstar._system.setVirtualSite(vid, vsite)

    def is_instance(self, filter):
        return filter in {"virtual_site", "vsite", "2fd"}
