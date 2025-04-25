from .force import Force
import openmm as mm  # type: ignore[import-untyped]


class HarmonicBond(Force):
    _members = 2

    def _set_force_obj(self):
        self._force_obj = mm.HarmonicBondForce()

    def _add_to_force_obj(self, params):
        self._force_obj.addBond(*params)

    def is_instance(self, filter):
        return filter in {"bond", "harmonic_bond"}
