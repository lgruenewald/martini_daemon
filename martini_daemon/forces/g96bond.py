from .force import Force
import openmm as mm


class G96Bond(Force):
    _members = 2

    def _set_force_obj(self):
        # Fourth power G96 potential
        self._force_obj = mm.CustomBondForce(
            "0.25 * k * (r^2-b^2)^2"
        )
        self._force_obj.addPerBondParameter("b")  # equilibrium length
        self._force_obj.addPerBondParameter("k")  # force constant

    def _add_to_force_obj(self, params):
        i, j, length, force = params
        self._force_obj.addBond(i, j, [length, force])

    def is_instance(self, filter):
        return filter in {"bond", "g96_bond"}
