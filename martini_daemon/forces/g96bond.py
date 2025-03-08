from .force import Force
import openmm as mm


class G96Bond(Force):
    def _build(self):
        # Fourth power G96 potential
        self._force_obj = mm.CustomBondForce(
            "0.25 * k * (r^2-b^2)^2"
        )
        self._force_obj.addPerBondParameter("b")  # equilibrium length
        self._force_obj.addPerBondParameter("k")  # force constant
        for (i, j, length, force) in filter(None, self._list):
            self._force_obj.addBond(i, j, [length, force])

    def add(self, members, params):
        i, j = members
        length, force = params
        self._list.append((i, j, length, force))
        if not self._rebuild:
            self._force_obj.addBond(i, j, [length, force])
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, *_ = self._list[id]
        return [i, j]

    def update_params(self, id, length, kb, kcub):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter in {"bond", "g96_bond"}
