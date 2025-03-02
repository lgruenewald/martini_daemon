import openmm as mm
from .force import Force


class QuarticAngle(Force):
    def _build(self):
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
        for (i, j, k, angle, c0, c1, c2, c3, c4) in filter(None, self._list):
            self._force_obj.addAngle(i, j, k, (angle, c0, c1, c2, c3, c4))

    def add(self, i, j, k, angle, c0, c1, c2, c3, c4):
        self._list.append((i, j, k, angle, c0, c1, c2, c3, c4))
        if not self._rebuild:
            self._force_obj.addAngle(i, j, k, (angle, c0, c1, c2, c3, c4))
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, *_ = self._list[id]
        return [i, j, k]

    def update_params(self, id, theta, force):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter in {"angle", "quartic_angle"}
