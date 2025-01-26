from .force import Force
import openmm as mm


class HarmonicBond(Force):
    """All bonds in the system
    indices are called bond_id
    values are (i: part_id, j: part_id, length: float, force: float)"""

    def _build(self):
        self._force_obj = mm.HarmonicBondForce()
        for (i, j, length, force) in filter(None, self._list):
            self._force_obj.addBond(i, j, length, force)

    def add(self, i, j, length, force):
        self._list.append((i, j, length, force))
        if not self._rebuild:
            self._force_obj.addBond(i, j, length, force)
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, _, _ = self._list[id]
        return [i, j]

    def update_params(self, id, length, force):
        raise NotImplementedError

