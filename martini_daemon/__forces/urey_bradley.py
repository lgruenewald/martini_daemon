import openmm as mm
from .force import Force


class UreyBradley(Force):
    _members = 3

    def _set_force_obj(self):
        self._force_obj = mm.CustomCompoundBondForce(
            3,  # 3 particles per compund bond force
            "0.5*k*((angle(p1,p2,p3)-theta0)^2)+"  # angle part
            "0.5*kUB*((distance(p1,p3)-r13)^2)"  # distance part
        )
        self._force_obj.addPerBondParameter("theta0")
        self._force_obj.addPerBondParameter("k")
        self._force_obj.addPerBondParameter("r13")
        self._force_obj.addPerBondParameter("kUB")

    def _add_to_force_obj(self, params):
        i, j, k, theta, force, r13, k_UB = params
        self._force_obj.addBond((i, j, k), (theta, force, r13, k_UB))

    _filters = {"angle", "urey_bradley"}
