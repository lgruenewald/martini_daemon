from abc import ABCMeta, abstractmethod
import openmm as mm
#from .system import System


class Force(metaclass=ABCMeta):

    def __init__(self, system):# System):
        self.force: mm.Force | None = None
        self.__force_group: int | None = None
        self.system = system
        # If rebuild is set to True, it means that self.force is  no longer valid, or does not exist.
        # Happens e.g. during construction or if a bond was removed.

    def build(self) -> None:
        if self.force is None and self.should_build():
            self.set_force_obj()
            assert self.force is not None
            self.prepare_force_obj()
            self.system.add_mm_force(self.force)

    def destroy(self) -> None:
        self.system.remove_mm_force(self.force)
        self.force = None

    def prepare_force_obj(self) -> None:
        self.force.setName(self.get_name())

    def should_build(self) -> bool:
        """
        May be used by bonded forces e.g. to indicate whether there is anything added to it.
        """
        return True

    @classmethod
    @abstractmethod
    def is_coupling(cls):
        raise NotImplementedError

    @abstractmethod
    def delta_degrees_of_freedom(self) -> int:
        raise NotImplementedError

    @abstractmethod
    def set_force_obj(self) -> None:
        """
        Define self.force, set it to an OpenMM force object. This will get added to the OpenMM System by the
        martini_daemon.System that owns this Force.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def get_name(cls) -> str:
        """
        There can only be one force with this name per martini_daemon.System.
        """
        raise NotImplementedError

    def has_force_obj(self) -> bool:
        return self.force is not None

    # TODO - reinitialize handling, atom type/charge change handling