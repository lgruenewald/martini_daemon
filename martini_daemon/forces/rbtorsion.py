from .force import Force
import openmm as mm


class RBTorsion(Force):
    def _build(self):
        self._force_obj = mm.RBTorsionForce()
        for (i, j, k, l, c0, c1, c2, c3, c4, c5) in filter(None, self._list):
            self._force_obj.addTorsion(i, j, k, l, c0, c1, c2, c3, c4, c5)

    def add(self, members, params):
        i, j, k, l = members
        c0, c1, c2, c3, c4, c5 = params
        self._list.append((i, j, k, l, c0, c1, c2, c3, c4, c5))
        if not self._rebuild:
            self._force_obj.addTorsion(i, j, k, l, c0, c1, c2, c3, c4, c5)
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, l, _, _, _ = self._list[id]
        return [i, j, k, l]

    def update_params(self, id, theta, force, mult):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter in {"dihedral", "rb_torsion"}
