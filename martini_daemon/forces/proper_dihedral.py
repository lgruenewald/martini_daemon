from .force import Force
import openmm as mm


class ProperDihedral(Force):
    _members = 4

    def _set_force_obj(self):
        self._force_obj = mm.PeriodicTorsionForce()

    def _add_to_force_obj(self, params):
        i, j, k, l, theta, force, mult = params
        self._force_obj.addTorsion(i, j, k, l, mult, theta, force)

    _filters = {"dihedral", "proper_dihedral"}
