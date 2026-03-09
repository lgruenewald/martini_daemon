from .force import Force
import openmm as mm


class CubicBond(Force):
    _members = 2

    def _set_force_obj(self):
        self._force_obj = mm.CustomBondForce(
            "kb * (r - b)^2 + kb * kcub * (r - b)^3"
        )
        self._force_obj.addPerBondParameter("b")  # equilibrium length
        self._force_obj.addPerBondParameter("kb")  # force constant
        self._force_obj.addPerBondParameter("kcub")  # cubic force constant

    def _add_to_force_obj(self, params):
        i, j, length, kb, kcub = params
        self._force_obj.addBond(i, j, [length, kb, kcub])

    _filters = {"bond", "cubic_bond"}
