from .force import Force
import openmm as mm
import math


class ImproperDihedral(Force):
    _members = 4

    def _set_force_obj(self):
        self._force_obj = mm.CustomTorsionForce(
            "0.5*k*(thetap-theta0)^2;"
            "thetap = step(-plus)*2*pi+theta+step(minus)*(-2*pi);"
            "plus=theta+pi-theta0;"
            "minus=theta-pi-theta0;"
            f"pi = {math.pi:.14f}"
        )
        self._force_obj.addPerTorsionParameter("theta0")
        self._force_obj.addPerTorsionParameter("k")

    def _add_to_force_obj(self, params):
        i, j, k, l, theta, force = params
        self._force_obj.addTorsion(i, j, k, l, (theta, force))

    _filters = {"dihedral", "improper_dihedral"}
