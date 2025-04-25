from ..forces.force import Interaction


class VirtualSite():
    # A parent class for virtual sites to use
    # Same interface as Force, so Interactions can be used

    def __init__(self, sysstar):
        self._list = []
        self._sysstar = sysstar
        self._built = False

    # overriden in child classes
    def _make_vsite(self, vid, members, params):
        raise NotImplementedError

    def is_instance(self, filter) -> bool:
        raise NotImplementedError

    def add(self, members, params) -> Interaction:
        if self._built:
            raise ValueError("Virtual sites cannot be added during the simulation.")
        vid, *other = members
        self._list.append((vid, other, params))
        self._sysstar.vsites.append(vid)
        return Interaction(self, len(self._list) - 1)

    def get_members(self, i):
        vid, members, _ = self._list[i]
        return [vid, *members]

    def remove(self, i):
        raise ValueError("Virtual sites cannot be removed during the simulation.")

    def build(self) -> None:
        if not self._built:
            self._built = True
            for vid, members, params in self._list:
                self._make_vsite(vid, members, params)

