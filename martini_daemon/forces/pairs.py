from .force import Force
import openmm as mm


class Pairs(Force):
    _members = 2

    types: dict[tuple[str, str], tuple[float, float]]

    def __init__(self, sysstar):
        super().__init__(sysstar)
        self.types = {}

    def add_type(self, type1, type2, c6, c12):
        self.types[(type1, type2)] = (c6, c12)
        self.types[(type2, type1)] = (c6, c12)

    def _set_force_obj(self):
        epsilon_r = self._sysstar.epsilon_r
        self._force_obj = mm.CustomBondForce(
            "LJ + ES;"
            "LJ = (C12 / r^12 - C6 / r^6);"
            "ES = f*qprod/epsilon_r/r;"
            f"epsilon_r = {epsilon_r};"
            "f = 138.935458;"
        )
        self._force_obj.addPerBondParameter("qprod")
        self._force_obj.addPerBondParameter("C6")
        self._force_obj.addPerBondParameter("C12")

    def _add_to_force_obj(self, params):
        if len(params) == 2:
            params.append(None)
            params.append(None)

        i, j, p1, p2 = params
        t1, q1, _ = self._sysstar.get_particle_details(i)
        t2, q2, _ = self._sysstar.get_particle_details(j)
        qprod = q1 * q2
        if p1 is None or p2 is None:
            p1, p2 = self.types[(t1, t2)]
        self._force_obj.addBond(i, j, [qprod, p1, p2])

    def is_instance(self, filter):
        return filter == "pair"
