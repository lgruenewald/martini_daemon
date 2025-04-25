import openmm as mm
from .force import Force


class CrossBondBond(Force):
    _members = 3

    def _set_force_obj(self):
        self._force_obj = mm.CustomCompoundBondForce(
            3,  # 3 particles per compund bond force
            "k*(distance(p1,p2)-r1)*(distance(p3,p2)-r2)"
        )
        self._force_obj.addPerBondParameter("r1")
        self._force_obj.addPerBondParameter("r2")
        self._force_obj.addPerBondParameter("k")

    def _add_to_force_obj(self, params):
        i, j, k, r1, r2, force = params
        self._force_obj.addBond((i, j, k), (r1, r2, force))

    def is_instance(self, filter):
        return filter in {"angle", "cross_bond_bond"}
