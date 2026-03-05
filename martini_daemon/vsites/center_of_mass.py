import openmm as mm
from .vsite import VirtualSite


class VSiteCenterOfMass(VirtualSite):
    def _make_vsite(self, vid, members, params):
        n = len(members)
        masses = []
        sum = 0.
        for i in members:
            _, _, m, _, _ = self._sysstar.get_atom_details(i)
            masses.append(m)
            sum += m
        weights = [m/sum for m in masses]
        vsite = mm.LocalCoordinatesSite(
            members,  # atoms
            weights,  # origin weights
            [0.0] * n,  # x direction weight
            [0.0] * n,  # y direction weight
            [0.0, 0.0, 0.0]  # coordinates
        )
        self._sysstar._system.setVirtualSite(vid, vsite)

    _filters = {"virtual_site", "vsite", "com"}
