from .force import Force
import openmm as mm


class RBTorsion(Force):
    _members = 4

    def _set_force_obj(self):
        self._force_obj = mm.RBTorsionForce()

    def _add_to_force_obj(self, params):
        i, j, k, l, c0, c1, c2, c3, c4, c5 = params
        self._force_obj.addTorsion(i, j, k, l, c0, c1, c2, c3, c4, c5)

    def is_instance(self, filter):
        return filter in {"dihedral", "rb_torsion"}
