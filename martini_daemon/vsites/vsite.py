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

    _filters = set()

    _destroyable = False

    def is_instance(self, filter) -> bool:
        return filter in self._filters

    def add(self, members, params) -> Interaction:
        if self._built:
            raise ValueError(
                "Virtual sites cannot be added during the simulation."
            )
        vid, *other = members
        self._list.append((vid, other, params))
        self._sysstar.vsites.append(vid)
        return Interaction(self, len(self._list) - 1)

    def get_members(self, i):
        vid, members, _ = self._list[i]
        return [vid, *members]

    def remove(self, i):
        raise ValueError(
            "Virtual sites cannot be removed during the simulation."
        )

    def build(self) -> None:
        if not self._built:
            self._built = True
            for vid, members, params in self._list:
                self._make_vsite(vid, members, params)

    def get_index(self) -> int | None:
        for i, f in enumerate(self._sysstar.modular_forces):
            if f == self:
                return i
        return None

    def get_class_name(self) -> str:
        return type(self).__name__

    def save(self, f) -> None:
        f.dump(self.get_class_name())
        f.dump(self.get_index())
        f.dump(self._list)

    def load(self, f) -> None:
        if self._built:
            raise ValueError("Can only load before _built")
        assert f.load() == self.get_class_name(), (
            "Attempt to load a checkpoint from a different version of daemon."
        )
        assert f.load() == self.get_index(), (
            "Attempt to load a checkpoint from a different version of daemon."
        )
        self._list = f.load()

    def __len__(self) -> int:
        return len(self._list)

    def set_force_group(self, fg) -> bool:
        return False

    def get_force_group(self) -> None:
        return None

