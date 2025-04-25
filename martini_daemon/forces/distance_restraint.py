from .force import Force
import openmm as mm


class DistanceRestraint(Force):
    _members = 2

    def _set_force_obj(self):
        self._force_obj = mm.CustomBondForce(
            "Vlow*step(low-r)+Vup1*step(r-up1)*step(up2-r)+Vup2*step(r-up2); "
            "Vlow = 0.5 * k *(r-low)^2; "
            "Vup1 = 0.5 * k * (r-up1)^2; "
            "Vup2 = 0.5 * k * (up2 - up1) * (2*r-up2-up1);"
        )
        self._force_obj.addPerBondParameter("low")
        self._force_obj.addPerBondParameter("up1")
        self._force_obj.addPerBondParameter("up2")
        self._force_obj.addPerBondParameter("k")

    def _add_to_force_obj(self, params):
        i, j, low, up1, up2, force = params
        self._force_obj.addBond(i, j, [low, up1, up2, force])

    def is_instance(self, filter):
        return filter in {"bond", "distance_restraint"}
