from __future__ import annotations
import openmm as mm
from dataclasses import dataclass


class Force():
    # if rebuild is set to True, it means that the force object _force_obj is
    # no longer valid, or does not exist. If rebuild is false, as much effort
    # should be done to keep _force_obj updated as possible
    # rebuild is True for example before force_obj is built in the first place,
    # or when removing interactions

    def __init__(self, sysstar):
        self._list = []
        self._sysstar = sysstar
        self._rebuild: bool = True
        self._force_obj: mm.Force = None

    # classes inheriting force must set these three functions and one int
    def _set_force_obj(self) -> None:
        raise NotImplementedError(f"_set_force_obj is not implemented for {self}")

    def _add_to_force_obj(self, params) -> None:
        raise NotImplementedError

    def is_instance(self, filter: str) -> bool:
        raise NotImplementedError

    _members: int

    def add(self, members: list[int] | tuple[int], params: list[float] | tuple[float]) -> Interaction:
        params = tuple(members) + tuple(params)
        self._list.append(params)
        if not self._rebuild:
            self._add_to_force_obj(params)
            self._sysstar._reinitialize = True
        return Interaction(self, len(self._list) - 1)

    def get_members(self, i: int) -> list[int]:
        params = self._list[i]
        return params[:self._members]

    def remove(self, i: int) -> None:
        self._list[i] = None
        self._rebuild = True

    def build(self) -> None:
        if self._rebuild:
            self.destroy()
            self._set_force_obj()
            for params in filter(None, self._list):
                self._add_to_force_obj(params)
            self._force_obj.setUsesPeriodicBoundaryConditions(True)
            self._rebuild = False
            self._sysstar._forces_list.append(self._force_obj)
            self._sysstar._system.addForce(self._force_obj)
            self._sysstar._reinitialize = True

    def destroy(self) -> bool:
        if self._force_obj is None:
            return False
        for i, f in enumerate(self._sysstar._forces_list):
            if f == self._force_obj:
                del self._sysstar._forces_list[i]
                self._sysstar._system.removeForce(i)
                self._sysstar._reinitialize = True
                return True
        return False


@dataclass
class Interaction():
    _force: Force
    _index: int

    def get_members(self) -> list[int]:
        return self._force.get_members(self._index)

    def remove(self) -> None:
        self._force.remove(self._index)

    def is_instance(self, filter: str) -> bool:
        return self._force.is_instance(filter)

    def __hash__(self):
        return hash((self._index, self._force.__class__))

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self._force == other._force and self._index == other._index
