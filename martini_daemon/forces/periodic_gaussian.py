from .force import Force
import openmm as mm
import math

# TODO split into plugin and own repo
class PeriodicGaussian(Force):
    _members = 4

    def _set_force_obj(self):
        self._force_obj = mm.CustomTorsionForce(
            "-A*exp(-k*x^2);"
            "x=2*pi*((theta-theta0)/(2*pi)-floor((theta-theta0+pi)/(2*pi)));"
            f"pi={math.pi};"
        )
        self._force_obj.addPerTorsionParameter("theta0")
        self._force_obj.addPerTorsionParameter("A")
        self._force_obj.addPerTorsionParameter("k")

    def _add_to_force_obj(self, params):
        i, j, k, L, theta, depth, force = params
        self._force_obj.addTorsion(i, j, k, L, [theta, depth, force])

    _filters = {"dihedral", "periodic_gaussian"}
