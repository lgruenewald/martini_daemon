from .force import Force
import openmm as mm


class RestrictedDihedral(Force):
    _members = 4

    def _set_force_obj(self):
        self._force_obj = mm.CustomTorsionForce(
            "0.5*k*((cos(theta)-cos(theta0))^2)/((sin(theta))^2)"
        )
        self._force_obj.addPerTorsionParameter("theta0")
        self._force_obj.addPerTorsionParameter("k")

    def _add_to_force_obj(self, params):
        i, j, k, l, theta, force = params
        self._force_obj.addTorsion(i, j, k, l, [theta, force])

    _filters = {"dihedral", "restricted_dihedral"}
