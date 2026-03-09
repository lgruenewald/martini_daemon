import openmm as mm
from .force import Force


class LinearAngle(Force):
    _members = 3

    def _set_force_obj(self):
        self._force_obj = mm.CustomCompoundBondForce(
            3,  # 3 particles per compund bond force
            "0.5*k*distj2; "
            "distj2=(xj-x2)^2+(yj-y2)^2+(zj-z2)^2; "
            "xj=a*x1+(1-a)*x3; "
            "yj=a*y1+(1-a)*y3; "
            "zj=a*z1+(1-a)*z3;"
        )
        self._force_obj.addPerBondParameter("a")
        self._force_obj.addPerBondParameter("k")

    def _add_to_force_obj(self, params):
        i, j, k, a, force = params
        self._force_obj.addBond((i, j, k), (a, force))

    _filters = {"angle", "linear_angle"}
