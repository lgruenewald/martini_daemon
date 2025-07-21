import openmm as mm
from .force import Force


class CombinedBendingTorsion(Force):
    _members = 4

    def _set_force_obj(self):
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

    def _add_to_force_obj(self, params):
        members = params[:4]
        params = params[4:]
        self._force_obj.addBond(members, params)

    _filters = {"dihedral", "combined_bending_torsion"}
