import openmm as mm
from ..forces.force import Force, Interaction


class VSite4fdn(Force):

    def _make_vsite(self, vid, i, j, k, l, a, b, c):
        vsite = mm.LocalCoordinatesSite(
            [i, j, k, l],  # particles
            [1.0, 0.0, 0.0, 0.0],  # origin weights
            [1.0 - a, -1.0, a, 0.0],  # x direction weight
            [1.0 - b, -1.0, 0.0, b],  # y direction weight
            [0.0, 0.0, c]  # coordinates
        )
        self._sysstar._system.setVirtualSite(vid, vsite)

    def _build(self):
        for (vid, i, j, k, l, a, b, c) in self._list:
            self._make_vsite(vid, i, j, k, l, a, b, c)

    def add(self, members, params):
        vid, i, j, k, l = members
        a, b, c = params
        self._list.append((vid, i, j, k, l, a, b, c))
        self._sysstar.vsites.append(vid)
        if not self._rebuild:
            raise Exception("Can't add virtual sites during run")
        return self._interaction()

    def get_members(self, i):
        vid, i, j, k, l, _,  _, _ = self._list[i]
        return [vid, i, j, k, l]

    def update_params(self, i, *params):
        raise NotImplementedError

    def remove(self, i):
        raise NotImplementedError

    def build(self):
        # only should get called once when building it initially
        if self._rebuild:
            self._build()
            self._rebuild = False
            self._sysstar._reinitialize = True

    def destroy(self):
        # no _force_obj, so it shold never be called like this
        raise NotImplementedError

    def _interaction(self):
        return Interaction(self, len(self._list) - 1)

