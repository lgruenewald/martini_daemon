import openmm as mm
from .force import Force


class UreyBradley(Force):
    def _build(self):
        self._force_obj = mm.CustomCompoundBondForce(
            3,  # 3 particles per compund bond force
            "0.5*k*((angle(p1,p2,p3)-theta0)^2)+"  # angle part
            "0.5*kUB*((distance(p1,p3)-r13)^2)"  # distance part
        )
        self._force_obj.addPerBondParameter("theta0")
        self._force_obj.addPerBondParameter("k")
        self._force_obj.addPerBondParameter("r13")
        self._force_obj.addPerBondParameter("kUB")
        for (i, j, k, theta, force, r13, k_UB) in filter(None, self._list):
            self._force_obj.addBond((i, j, k), (theta, force, r13, k_UB))

    def add(self, i, j, k, theta, force, r13, k_UB):
        self._list.append((i, j, k, theta, force, r13, k_UB))
        if not self._rebuild:
            self._force_obj.addBond((i, j, k), (theta, force, r13, k_UB))
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, *_ = self._list[id]
        return [i, j, k]

    def update_params(self, id, theta, force):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter in {"angle", "urey_bradley"}
