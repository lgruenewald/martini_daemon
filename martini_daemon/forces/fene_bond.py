from .force import Force
import openmm as mm


class FENEBond(Force):
    """All fene (finitely extensible nonlinear elastic) bonds in the system
    indices are called bond_id
    values are (i: part_id, j: part_id, length, D, beta: float)"""

    def _build(self):
        self._force_obj = mm.CustomBondForce(
            "- 0.5 * k * b^2 * log(1 - r^2 / b^2)"
        )
        self._force_obj.addPerBondParameter("b")  # equilibrium length
        self._force_obj.addPerBondParameter("k")  # force constant
        for (i, j, length, k) in filter(None, self._list):
            self._force_obj.addBond(i, j, [length, k])

    def add(self, members, params):
        i, j = members
        length, k = params
        self._list.append((i, j, length, k))
        if not self._rebuild:
            self._force_obj.addBond(i, j, [length, k])
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, _, _ = self._list[id]
        return [i, j]

    def update_params(self, id, length, kb, kcub):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter in {"bond", "fene_bond"}

