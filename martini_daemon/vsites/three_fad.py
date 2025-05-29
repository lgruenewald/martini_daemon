import openmm as mm
from .vsite import VirtualSite
from math import cos, sin


class VSite3fad(VirtualSite):

    def _make_vsite(self, vid, members, params):
        i, j, k = members
        theta, d = params
        vsite = mm.LocalCoordinatesSite(
            [i, j, k],  # atoms
            [1.0, 0.0, 0.0],  # origin weights
            [-0.5, 0.5, 0.0],  # x direction weight
            [0.0, -0.5, 0.5],  # y direction weight
            [d * cos(theta), d * sin(theta), 0.0]  # coordinates
        )
        self._sysstar._system.setVirtualSite(vid, vsite)

    def is_instance(self, filter):
        return filter in {"virtual_site", "vsite", "3fad"}
