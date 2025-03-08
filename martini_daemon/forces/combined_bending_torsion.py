import openmm as mm
from .force import Force


class CombinedBendingTorsion(Force):
    def _build(self):
        self._force_obj = mm.CustomCompoundBondForce(
            4,
            "k*sintheta0^3*sintheta1^3*(a0 + a1*cosphi + a2*cosphi^2 + a3*cosphi^3 + a4*cosphi^4); "
            "sintheta0 = sin(angle(p1, p2, p3));"
            "sintheta1 = sin(angle(p2, p3, p4));"
            "cosphi = cos(dihedral(p1, p2, p3, p4));",
        )
        self._force_obj.addPerBondParameter("k")
        self._force_obj.addPerBondParameter("a0")
        self._force_obj.addPerBondParameter("a1")
        self._force_obj.addPerBondParameter("a2")
        self._force_obj.addPerBondParameter("a3")
        self._force_obj.addPerBondParameter("a4")
        for (i, j, k, l, force, a0, a1, a2, a3, a4) in filter(None, self._list):
            self._force_obj.addBond((i, j, k, l), (force, a0, a1, a2, a3, a4))

    def add(self, members, params):
        i, j, k, l = members
        force, a0, a1, a2, a3, a4 = params
        self._list.append((i, j, k, l, force, a0, a1, a2, a3, a4))
        if not self._rebuild:
            self._force_obj.addBond((i, j, k, l), (force, a0, a1, a2, a3, a4))
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, k, l, _, _, _, _, _, _ = self._list[id]
        return [i, j, k, l]

    def update_params(self, id, theta, force):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter in {"dihedral", "combined_bending_torsion"}
