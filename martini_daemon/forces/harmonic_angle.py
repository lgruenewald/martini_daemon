from .force import Force
import openmm as mm

class HarmonicAngle(Force):
    """All the angles in the system
    indices are called angle_id
    values are (i: part_id, j: part_id, k: part_id, theta: float, force: float)
    """

    def _build(self):
        self._force_obj = mm.HarmonicAngleForce()
        for (i, j, k, theta, force) in filter(None, self._list):
            self._force_obj.addAngle(i, j, k, theta, force)

    def add(self, i, j, k, theta, force):
        self._list.append((i, j, k, theta, force))
        if not self._rebuild:
            self._force_obj.addAngle(i, j, k, theta, force)
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, _, _ = self._list[id]
        return [i, j, k]

    def update_params(self, id, theta, force):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter in {"angle", "harmonic_angle"}
