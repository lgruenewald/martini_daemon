import openmm as mm  # type: ignore[import-untyped]
from ..forces.force import Force, Interaction


class VSiteCenterOfMass(Force):
    # TODO document limitation about updating of masses during reaction
    def _make_vsite(self, vid, members):
        n = len(members)
        masses = []
        sum = 0.
        for i in members:
            _, _, m = self._sysstar.get_particle_details(i)
            masses.append(m)
            sum += m
        weights = [m/sum for m in masses]
        vsite = mm.LocalCoordinatesSite(
            members,  # particles
            weights,  # origin weights
            [0.0] * n,  # x direction weight
            [0.0] * n,  # y direction weight
            [0.0, 0.0, 0.0]  # coordinates
        )
        self._sysstar._system.setVirtualSite(vid, vsite)

    def _build(self):
        for (vid, members) in self._list:
            self._make_vsite(vid, members)

    def add(self, members, params):
        # params not used
        vid, *members = members
        self._list.append((vid, members))
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

    def is_instance(self, filter):
        return filter in {"virtual_site", "vsite", "com"}

