from .force import Force
import openmm as mm


class MorseBond(Force):
    _members = 2

    def _set_force_obj(self):
        self._force_obj = mm.CustomBondForce(
            "D * (1 - exp(-beta * (r - b)))^2"
        )
        self._force_obj.addPerBondParameter("b")  # equilibrium length
        self._force_obj.addPerBondParameter("D")  # force constant
        self._force_obj.addPerBondParameter("beta")  # cubic force constant

    def _add_to_force_obj(self, params):
        i, j, length, D, beta = params
        self._force_obj.addBond(i, j, [length, D, beta])

    _filters = {"bond", "morse_bond"}
