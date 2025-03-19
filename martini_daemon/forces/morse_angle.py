import openmm as mm
from .force import Force


class MorseAngle(Force):
    def _build(self):
        self._force_obj = mm.CustomCompoundBondForce(
            3,  # 3 particles per compund bond force
            "D*(1-morse)^2+"  # bond part
            "k*dscale*(angle(p1, p2, p3)-theta0)^2; "  # angle part
            "dscale=2*morse/(1+morse); "
            "morse=exp(-beta*(distance(p1, p2) - r0));"
        )
        self._force_obj.addPerBondParameter("theta0")
        self._force_obj.addPerBondParameter("k")
        self._force_obj.addPerBondParameter("r0")
        self._force_obj.addPerBondParameter("D")
        self._force_obj.addPerBondParameter("beta")
        for (i, j, k, theta, force, r0, D, beta) in filter(None, self._list):
            self._force_obj.addBond((i, j, k), (theta, force, r0, D, beta))

    def add(self, members, params):
        i, j, k = members
        theta, force, r0, D, beta = params
        self._list.append((i, j, k, theta, force, r0, D, beta))
        if not self._rebuild:
            self._force_obj.addBond((i, j, k), (theta, force, r0, D, beta))
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, *_ = self._list[id]
        return [i, j, k]

    def update_params(self, id, theta, force):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter in {"angle", "morse_angle"}
