from .force import Force
import openmm as mm


class HarmonicAngle(Force):
    _members = 3

    def _set_force_obj(self):
        self._force_obj = mm.HarmonicAngleForce()

    def _add_to_force_obj(self, params):
        self._force_obj.addAngle(*params)

    _filters = {"angle", "harmonic_angle"}
