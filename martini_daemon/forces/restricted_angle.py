import openmm as mm
from .force import Force


class RestrictedAngle(Force):
    def _build(self):
        self._force_obj = mm.CustomAngleForce(
            "0.5*k*(cos(theta)-cos(theta0))^2/sin(theta)^2"
        )
        self._force_obj.addPerAngleParameter("theta0")
        self._force_obj.addPerAngleParameter("k")
        for (i, j, k, theta, force) in filter(None, self._list):
            self._force_obj.addAngle(i, j, k, (theta, force))

    def add(self, members, params):
        i, j, k = members
        theta, force = params
        self._list.append((i, j, k, theta, force))
        if not self._rebuild:
            self._force_obj.addAngle(i, j, k, (theta, force))
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, _, _ = self._list[id]
        return [i, j, k]

    def update_params(self, id, theta, force):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter in {"angle", "restricted_angle"}
