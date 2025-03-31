from .force import Force
import openmm as mm  # type: ignore[import-untyped]


class ProperDihedral(Force):
    """All the proper dihedrals in the system
    indices are called proper_dihedral_id
    values are (i: part_id, j: part_id, k: part_id, l: part_id,
     theta: float, force: float, multiplicity: int)"""

    def _build(self):
        self._force_obj = mm.PeriodicTorsionForce()
        for (i, j, k, l, theta, force, mult) in filter(None, self._list):
            self._force_obj.addTorsion(i, j, k, l, mult, theta, force)

    def add(self, members, params):
        i, j, k, l = members
        theta, force, mult = params
        self._list.append((i, j, k, l, theta, force, mult))
        if not self._rebuild:
            self._force_obj.addTorsion(i, j, k, l, mult, theta, force)
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, l, _, _, _ = self._list[id]
        return [i, j, k, l]

    def update_params(self, id, theta, force, mult):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter in {"dihedral", "proper_dihedral"}
