from .force import Force
import openmm as mm
import math


class ImproperDihedral(Force):
    """All the improper dihedrals in the system
    indices are called improper_dihedral_id
    values are (i: part_id, j: part_id, k: part_id, l: part_id,
     theta: float, force: float)
    """

    def _build(self):
        self._force_obj = mm.CustomTorsionForce(
            "0.5*k*(thetap-theta0)^2;"
            "thetap = step(-plus)*2*pi+theta+step(minus)*(-2*pi);"
            "plus=theta+pi-theta0;"
            "minus=theta-pi-theta0;"
            f"pi = {math.pi:.14f}"
        )
        self._force_obj.addPerTorsionParameter("theta0")
        self._force_obj.addPerTorsionParameter("k")
        for (i, j, k, l, theta, force) in filter(None, self._list):
            self._force_obj.addTorsion(i, j, k, l, (theta, force))

    def add(self, i, j, k, l, theta, force):
        self._list.append((i, j, k, l, theta, force))
        if not self._rebuild:
            self._force_obj.addTorsion(i, j, k, l, (theta, force))
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, l, _, _ = self._list[id]
        return [i, j, k, l]

    def update_params(self, id, theta, force):
        raise NotImplementedError

