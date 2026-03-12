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

    def build(self, must=False) -> None:
        """
        Call this before initializing or reinitializing the context. It actually creates the OpenMM force
        and adds it to the system.

        If must is True, then assume that the force object has been deleted and a new one is needed, don't check.
        Currently must=True is used for constraints <=> virtual sites
        """
        if must or self.should_build():
            self._set_force_obj()
            assert self.force is not None
            self._prepare_force_obj()
            self.system._add_mm_force(self.force)

    def _destroy(self) -> None:
        """
        Will destroy the OpenMM force and remove it from the OpenMM system.

        Note: the force will automatically get rebuilt before continuing the simulation.

        Protected because it should only be called by this class.
        """
        self.system._remove_mm_force(self.force)
        self.force = None

    def _prepare_force_obj(self) -> None:
        self.force.setName(self.get_name())

    def should_build(self) -> bool:
        """
        May be used by bonded forces e.g. to indicate whether there is anything added to it.
        """
        return self.force is None

    @classmethod
    @abstractmethod
    def is_coupling(cls):
        """
        Whether this Force wraps an OpenMM Force that is a type of coupling (temperature, pressure, COMM removal).
        """
        raise NotImplementedError

    @abstractmethod
    def delta_degrees_of_freedom(self) -> int:
        """
        How many degrees of freedom does this force remove from the system.

        For example, COMM removal should return -3.
        Virtual sites should return -3N where N is the number of virtual sites in the force.
        """
        raise NotImplementedError

    @abstractmethod
    def _set_force_obj(self) -> None:
        """
        Define self.force, set it to an OpenMM force object. This will get added to the OpenMM System by the
        martini_daemon.System that owns this Force.

        Protected because it should only be called by this class.
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

    def flag_atom_change(self, atom_id: int, change_charge: bool = False) -> None:
        pass