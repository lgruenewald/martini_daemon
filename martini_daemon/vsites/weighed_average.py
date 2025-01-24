import openmm as mm
from ..forces.force import Force, Interaction


class VSiteWeighedAverage(Force):
    # Force because it implements the same interface, but it reimplements
    # everything, the default helpers don't suit it well (no _force_obj).

    def _make_vsite(self, vid, members, weights):
        n = len(members)
        if n == 1:
            vsite = mm.LocalCoordinatesSite(
                [members[0], 0],
                [1.0, 0.0],
                [0.0, 0.0],
                [0.0, 0.0],
                [0.0, 0.0, 0.0]
            )
        elif n == 2:
            vsite = mm.TwoParticleAverageSite(
                members[0], members[1], weights[0], weights[1]
            )
        elif n == 3:
            vsite = mm.ThreeParticleAverageSite(
                members[0], members[1], members[2],
                weights[0], weights[1], weights[2]
            )
        else:
            vsite = mm.LocalCoordinatesSite(
                members,  # particles
                weights,  # origin weights
                [0.0] * n,  # x direction weight
                [0.0] * n,  # y direction weight
                [0.0, 0.0, 0.0]  # coordinates
            )
        self._sysstar._system.setVirtualSite(vid, vsite)

    def _build(self):
        for (vid, members, weights) in self._list:
            self._make_vsite(vid, members, weights)

    def add(self, members, weights):
        vid, *members = members
        self._list.append((vid, members, weights))
        self._sysstar.vsites.append(vid)
        if not self._rebuild:
            raise Exception("Can't add virtual sites during run")
        return self._interaction()

    def get_members(self, i):
        vid, members, _ = self._list[i]
        return [vid, *members]

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

