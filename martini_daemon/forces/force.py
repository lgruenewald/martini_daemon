import openmm as mm
from dataclasses import dataclass


class Force():
    _force_obj: mm.Force
    _list: list
    _rebuild: bool
    # if rebuild is set to True, it means that the force object _force_obj is
    # no longer valid, or does not exist. If rebuild is false, as much effort
    # should be done to keep _force_obj updated as possible
    # rebuild is True for example before force_obj is built in the first place,
    # or when removing elements from the bond

    def __init__(self, sysstar):
        self._list = []
        self._sysstar = sysstar
        self._rebuild = True
        self._force_obj = None

    def _build(self):
        pass

    def add(self, *params):
        pass

    def get_members(self, i):
        pass

    def update_params(self, i, *params):
        pass

    def remove(self, i):
        self._list[i] = None
        self._rebuild = True

    def build(self):
        if self._rebuild:
            self.destroy()
            self._build()
            self._force_obj.setUsesPeriodicBoundaryConditions(True)
            self._rebuild = False
            self._sysstar._forces_list.append(self._force_obj)
            self._sysstar._system.addForce(self._force_obj)
            self._sysstar._reinitialize = True

    def destroy(self):
        if self._force_obj is None:
            return False
        for i, f in enumerate(self._sysstar._forces_list):
            if f == self._force_obj:
                del self._sysstar._forces_list[i]
                self._sysstar._system.removeForce(i)
                self._sysstar._reinitialize = True
                return True
        return False

    def _interaction(self):
        return Interaction(self, len(self._list) - 1)


@dataclass
class Interaction():
    _force: Force
    _index: int

    def get_members(self):
        return self._force.get_members(self._index)

    def update_params(self, *params):
        return self._force.update_params(self._index, params)

    def remove(self):
        self._force.remove(self._index)

    def __hash__(self):
        return hash((self._index, self._force.__class__))

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self._force == other._force and self._index == other._index
