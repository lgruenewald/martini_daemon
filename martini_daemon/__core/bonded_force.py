from __future__ import annotations
import openmm as mm
from abc import ABC, abstractmethod
from .system import System

class BondedForce(ABC):

    def __init__(self, system: System):
        self.force = None
        self.force_group: int | None = None
        self.entries: dict[int, tuple[list[int], list[float]]] = {}
        self.system = system
        # if rebuild is set to True, it means that the force object _force_obj is
        # no longer valid, or does not exist. If rebuild is false, as much effort
        # should be done to keep _force_obj updated as possible
        # rebuild is True for example before force_obj is built in the first place,
        # or when removing interactions
        self.rebuild: bool = True

    @abstractmethod
    def set_force_obj(self) -> None:
        """
        Define self.force, set it to an OpenMM force object.
        """
        raise NotImplementedError

    def prepare_force_obj(self) -> None:
        """
        Prepare force object. By default, this sets the PBC handling to true. Override if this is not desired.
        """
        self.force.setUsesPeriodicBoundaryConditions(True)

    @abstractmethod
    def parse(self, members: list[int], params: list[float]) -> list[float]:
        """
        Given a list of members and parameters (as unwrapped = minimal pre-parsing), pa
        """
        raise NotImplementedError

    @abstractmethod
    def passes_filter(self, filter_: str) -> bool:
        return filter_ in self.filters

    @classmethod
    @abstractmethod
    def get_name(cls) -> str:
        raise NotImplementedError

    members: int
    filters: set[str] = set()

    def add(self, members: list[int] | tuple[int], params: list[float] | tuple[float]):# -> Interaction:
        params = tuple(members) + tuple(params)
        self._list.append(params)
        if not self._rebuild and self.has_force_obj():
            self._add_to_force_obj(params)
            self._sysstar._reinitialize = True
        #return Interaction(self, len(self._list) - 1)

    def get_members(self, i: int) -> None | list[int]:
        params = self._list[i]
        if params is None:
            return None
        return params[:self._members]

    def remove(self, i: int) -> None:
        self._list[i] = None
        self._rebuild = True

    def build(self) -> None:
        if self._rebuild:
            self.destroy()
            self._rebuild = False
            self._sysstar._reinitialize = True
            if len(self._list) > 0:
                self._set_force_obj()
                if self._force_group is not None:
                    self._force_obj.setForceGroup(self._force_group)
                for params in filter(None, self._list):
                    self._add_to_force_obj(params)
                if self._pbc:
                    self._force_obj.setUsesPeriodicBoundaryConditions(True)
                # TODO why not sysstar.add_force
                self._sysstar._forces_list.append(self._force_obj)
                self._sysstar._system.addForce(self._force_obj)

    def destroy(self) -> bool:
        if self._force_obj is None:
            return False
        for i, f in enumerate(self._sysstar._forces_list):
            if f == self._force_obj:
                self._force_obj = None
                del self._sysstar._forces_list[i]
                self._sysstar._system.removeForce(i)
                self._sysstar._reinitialize = True
                return True
        return False

    def get_index(self) -> int | None:
        for i, f in enumerate(self._sysstar.modular_forces):
            if f == self:
                return i
        return None


    def __len__(self) -> int:
        return len(self._list)

    def set_force_group(self, fg: int) -> bool:
        self._force_group = fg
        return True

    def get_force_group(self) -> int | None:
        return self._force_group

    def has_force_obj(self) -> bool:
        return self._force_obj is not None

