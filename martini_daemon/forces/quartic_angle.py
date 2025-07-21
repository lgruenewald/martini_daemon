import openmm as mm
from .force import Force


class QuarticAngle(Force):
    _members = 3

    def _set_force_obj(self):
        self._force_obj = mm.CustomAngleForce(
            "c0+"
            "c1*(theta-theta0)+"
            "c2*(theta-theta0)^2+"
            "c3*(theta-theta0)^3+"
            "c4*(theta-theta0)^4"
        )
        self._force_obj.addPerAngleParameter("theta0")
        self._force_obj.addPerAngleParameter("c0")
        self._force_obj.addPerAngleParameter("c1")
        self._force_obj.addPerAngleParameter("c2")
        self._force_obj.addPerAngleParameter("c3")
        self._force_obj.addPerAngleParameter("c4")

    def _add_to_force_obj(self, params):
        i, j, k, angle, c0, c1, c2, c3, c4 = params
        self._force_obj.addAngle(i, j, k, (angle, c0, c1, c2, c3, c4))

    _filters = {"angle", "quartic_angle"}
