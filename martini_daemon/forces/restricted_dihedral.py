from .force import Force
import openmm as mm


class RestrictedDihedral(Force):

    def _build(self):
        self._force_obj = mm.CustomTorsionForce(
            "0.5*k*((cos(theta)-cos(theta0))^2)/((sin(theta))^2)"
        )
        self._force_obj.addPerTorsionParameter("theta0")
        self._force_obj.addPerTorsionParameter("k")
        for (i, j, k, l, theta, force) in filter(None, self._list):
            self._force_obj.addTorsion(i, j, k, l, [theta, force])

    def add(self, i, j, k, l, theta, force):
        self._list.append((i, j, k, l, theta, force))
        if not self._rebuild:
            self._force_obj.addTorsion(i, j, k, l, [theta, force])
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, l, _, _, _ = self._list[id]
        return [i, j, k, l]

    def update_params(self, id, theta, force, mult):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter in {"dihedral", "restricted_dihedral"}
