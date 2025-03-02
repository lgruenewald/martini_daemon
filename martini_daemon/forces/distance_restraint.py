from .force import Force
import openmm as mm


class DistanceRestraint(Force):
    def _build(self):
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
        for (i, j, low, up1, up2, force) in filter(None, self._list):
            self._force_obj.addBond(i, j, [low, up1, up2, force])

    def add(self, i, j, low, up1, up2, force):
        self._list.append((i, j, low, up1, up2, force))
        if not self._rebuild:
            self._force_obj.addBond(i, j, [low, up1, up2, force])
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, *_ = self._list[id]
        return [i, j]

    def update_params(self, id, low, up1, up2, force):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter in {"bond", "distance_restraint"}
