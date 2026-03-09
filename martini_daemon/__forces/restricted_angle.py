import openmm as mm
from .force import Force


class RestrictedAngle(Force):
    _members = 3

    def _set_force_obj(self):
        self._force_obj = mm.CustomAngleForce(
            "0.5*k*(cos(theta)-cos(theta0))^2/sin(theta)^2"
        )
        self._force_obj.addPerAngleParameter("theta0")
        self._force_obj.addPerAngleParameter("k")

    def _add_to_force_obj(self, params):
        i, j, k, theta, force = params
        self._force_obj.addAngle(i, j, k, (theta, force))

    _filters = {"angle", "restricted_angle"}
