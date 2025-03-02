from .force import Force
import openmm as mm


class MorseBond(Force):
    """All morse bonds in the system
    indices are called bond_id
    values are (i: part_id, j: part_id, length, D, beta: float)"""

    def _build(self):
        self._force_obj = mm.CustomBondForce(
            "D * (1 - exp(-beta * (r - b)))^2"
        )
        self._force_obj.addPerBondParameter("b")  # equilibrium length
        self._force_obj.addPerBondParameter("D")  # force constant
        self._force_obj.addPerBondParameter("beta")  # cubic force constant
        for (i, j, length, D, beta) in filter(None, self._list):
            self._force_obj.addBond(i, j, [length, D, beta])

    def add(self, i, j, length, D, beta):
        self._list.append((i, j, length, D, beta))
        if not self._rebuild:
            self._force_obj.addBond(i, j, [length, D, beta])
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, *_ = self._list[id]
        return [i, j]

    def update_params(self, id, length, kb, kcub):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter in {"bond", "morse_bond"}
