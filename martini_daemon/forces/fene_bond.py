from .force import Force
import openmm as mm


class FENEBond(Force):
    _members = 2

    def _set_force_obj(self):
        self._force_obj = mm.CustomBondForce(
            "- 0.5 * k * b^2 * log(1 - r^2 / b^2)"
        )
        self._force_obj.addPerBondParameter("b")  # equilibrium length
        self._force_obj.addPerBondParameter("k")  # force constant

    def _add_to_force_obj(self, params):
        i, j, length, k = params
        self._force_obj.addBond(i, j, [length, k])

    def is_instance(self, filter):
        return filter in {"bond", "fene_bond"}

