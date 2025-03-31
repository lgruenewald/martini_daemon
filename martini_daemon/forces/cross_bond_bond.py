import openmm as mm  # type: ignore[import-untyped]
from .force import Force


class CrossBondBond(Force):
    def _build(self):
        self._force_obj = mm.CustomCompoundBondForce(
            3,  # 3 particles per compund bond force
            "k*(distance(p1,p2)-r1)*(distance(p3,p2)-r2)"
        )
        self._force_obj.addPerBondParameter("r1")
        self._force_obj.addPerBondParameter("r2")
        self._force_obj.addPerBondParameter("k")
        for (i, j, k, r1, r2, force) in filter(None, self._list):
            self._force_obj.addBond((i, j, k), (r1, r2, force))

    def add(self, members, params):
        i, j, k = members
        r1, r2, force = params
        self._list.append((i, j, k, r1, r2, force))
        if not self._rebuild:
            self._force_obj.addBond((i, j, k), (r1, r2, force))
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, *_ = self._list[id]
        return [i, j, k]

    def update_params(self, id, theta, force):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter in {"angle", "cross_bond_bond"}
