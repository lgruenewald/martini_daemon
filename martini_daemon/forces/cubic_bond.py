from .force import Force
import openmm as mm


class CubicBond(Force):
    """All cubic bonds in the system
    indices are called bond_id
    values are (i: part_id, j: part_id, length, kb, kcub: float)"""

    def _build(self):
        self._force_obj = mm.CustomBondForce(
            "kb * (r - b)^2 + kb * kcub * (r - b)^3"
        )
        self._force_obj.addPerBondParameter("b")  # equilibrium length
        self._force_obj.addPerBondParameter("kb")  # force constant
        self._force_obj.addPerBondParameter("kcub")  # cubic force constant
        for (i, j, length, kb, kcub) in filter(None, self._list):
            self._force_obj.addBond(i, j, [length, kb, kcub])

    def add(self, members, params):
        i, j = members
        length, kb, kcub = params
        self._list.append((i, j, length, kb, kcub))
        if not self._rebuild:
            self._force_obj.addBond(i, j, [length, kb, kcub])
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, _, _ = self._list[id]
        return [i, j]

    def update_params(self, id, length, kb, kcub):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter in {"bond", "cubic_bond"}

