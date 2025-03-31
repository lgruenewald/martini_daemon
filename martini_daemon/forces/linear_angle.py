import openmm as mm  # type: ignore[import-untyped]
from .force import Force


class LinearAngle(Force):
    def _build(self):
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
        for (i, j, k, a, force) in filter(None, self._list):
            self._force_obj.addBond((i, j, k), (a, force))

    def add(self, members, params):
        i, j, k = members
        a, force = params
        self._list.append((i, j, k, a, force))
        if not self._rebuild:
            self._force_obj.addBond((i, j, k), (a, force))
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, *_ = self._list[id]
        return [i, j, k]

    def update_params(self, id, theta, force):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter in {"angle", "linear_angle"}
