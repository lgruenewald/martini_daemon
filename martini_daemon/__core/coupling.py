# OpenMM coupling forces will be wrapped in this
# note: this class doesn't allow them to be removed later, as the settings to generate them are unknown


import openmm as mm

from .force import Force


def wrap_coupling(mm_force: mm.Force) -> type[Force]:
    """Given an OpenMM force that should act as a coupling for a simulation,
    it creates a Martini Daemon Force from it.
    """
    used = False

    class Coupling(Force):
        @classmethod
        def is_coupling(cls) -> bool:
            return True

        def delta_degrees_of_freedom(self) -> int:
            return 3 if isinstance(self._force, mm.CMMotionRemover) else 0

        def _set_force_obj(self) -> mm.Force:
            nonlocal used, mm_force
            assert not used, "wrapped couplings can be only set once"
            used = True

            return mm_force

        def _destroy(self) -> None:
            assert False, "wrapped couplings cannot be destroyed or reconstructed."

        @classmethod
        def get_name(cls) -> str:
            return type(mm_force).__name__.lower()

    return Coupling
