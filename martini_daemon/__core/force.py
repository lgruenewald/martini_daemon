from abc import ABCMeta, abstractmethod

import openmm as mm


class Force(metaclass=ABCMeta):
    # avoids typing to avoid circular dependency
    def __init__(self, system) -> None:  # noqa: ANN001
        """
        Parent metaclass for all forces in the system.

        Note the following invariant:

        - self._force should be not None if and only if it represents the most up-to-date state
            of the system and is present in the OpenMM system.
            When it becomes outdated it should be removed from System and self.force
            should be set to None.
        - if self._force is None, self._set_force_obj() will be called if
            self.should_build() returns True, before the energies/forces are read out.
            self._set_force_obj() is called by build(), which is in turn called
            by System automatically.

        You should implement the following in child classes:

        - get_name()
        - delta_degrees_of_freedom()
        - is_coupling()
        - _set_force_obj()

        You can optionally implement the following in child classes (if you know what you're doing):

        - _prepare_force_obj()
        - should_build()
        - flag_atom_change()
        - flag_add_atom()

        :param system: A Martini Daemon System instance. Note: type hint incorrect to avoid circular dependencies.
            self.system is set to system. Classes that inherit this class can safely import System, and can safely
            assume that self.system is always a Martini Daemon system.
        """
        self._force: mm.Force | None = None
        self.system = system

    def build(self, must: bool = False) -> None:
        """
        Build the OpenMM force associated with this Martini Daemon force.

        System calls this before initializing or reinitializing the context.
        It actually creates the OpenMM force and adds it to the system.

        :param must: If must is True, then assume that the force object has been deleted manually and a new one is
            needed, don't check.
            Currently must=True is only used for constraints <=> harmonic bonds conversion during minimization.
            Don't call with must=True if you don't know what you're doing.
        """
        if must or self.should_build():
            self._force = self._set_force_obj()
            if self._force is not None:
                self._prepare_force_obj()
                self.system._add_mm_force(self._force)

    def _destroy(self) -> None:
        """
        Will destroy the OpenMM force and remove it from the OpenMM system.

        Note: the force will automatically get rebuilt before continuing the simulation.

        Protected because it should only be called by methods in this class and its children.
        """
        if self._force is not None:
            self.system._remove_mm_force(self._force)
            self._force = None

    def _prepare_force_obj(self) -> None:
        """
        Prepare the OpenMM force object before it is added to System, after _set_force_obj().

        The default impl sets the name to the result of self.get_name(),
        which should be done, as System may rely on this.
        """
        assert self._force is not None
        self._force.setName(self.get_name())

    def should_build(self) -> bool:
        """
        Return whether the force should be built when appropriate.

        If it returns False, the force is either built or it does not need to be
        built (because it makes 0 difference whether it's present). By default
        returns True if self.force is None.

        May be used by bonded forces e.g. to indicate whether there is anything added to it.
        """
        return self._force is None

    @classmethod
    @abstractmethod
    def is_coupling(cls) -> bool:
        """Whether this Force wraps an OpenMM Force that is a type of coupling (temperature, pressure, COMM removal)."""
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
    def _set_force_obj(self) -> mm.Force:
        """
        Create the OpenMM force object.

        Define an OpenMM force object and return it. This will get added to the OpenMM System by the
        martini_daemon.System that owns this Force.

        This is called lazily, and may be called multiple times, therefore it should not have any other side effects
        than returning a force. It may read things from self.system though.

        Protected because it should only be called by this class.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def get_name(cls) -> str:
        """
        Return the name associated with this force.

        There can only be one force with this name per martini_daemon.System.
        """
        raise NotImplementedError

    def flag_atom_change(self, atom_id: int, change_charge: bool = False) -> None:
        """
        Let the force know that an atom type, charge or mass are changed.

        Called by System.
        Override if you need to update the OpenMM Force when this happens.

        Note: you should check if self.force exists first, before doing anything to it.
        """

    def flag_atom_add(self) -> None:
        """
        Let the force object know, that a new atom was added to the system.

        Called by System, but child classes can override.

        Note: currently this only happens during system construction.

        Note: if you inherit this, you should check if self.force exists, before doing anything to it.
        """

    def flag_add_exclusion(self, i: int, j: int) -> None:
        """
        Let the force object know, that a new exclusion was added to the system.

        Called by system, but child classes can override.

        If implementing a non bonded force, all nonbonded forces must have
        the same set of exclusions globally, due to OpenMM limitations.
        """

    def flag_remove_exclusion(self) -> None:
        """
        Let the force object know, that an exclusion was removed from the system.

        See flag_add_exclusion for further notes.
        """
