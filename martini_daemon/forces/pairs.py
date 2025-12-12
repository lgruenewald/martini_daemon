from .force import Force
import openmm as mm


class Pairs(Force):
    _members = 2

    types: dict[tuple[str, str], tuple[float, float]]

    def __init__(self, sysstar):
        super().__init__(sysstar)
        self.types = {}

    def add_type(self, type1, type2, sigma, epsilon):
        self.types[(type1, type2)] = (sigma, epsilon)
        self.types[(type2, type1)] = (sigma, epsilon)

    def _set_force_obj(self):
        epsilon_r = self._sysstar.nonbonded_force.epsilon_r
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
            params = (*params, None, None)

        i, j, sigma, epsilon = params
        t1, q1, _, _, _ = self._sysstar.get_atom_details(i)
        t2, q2, _, _, _ = self._sysstar.get_atom_details(j)
        qprod = q1 * q2
        if sigma is None or epsilon is None:
            if self.types.get((t1, t2)) is None:
                raise ValueError(
                    f"Unknown pair type ({t1}, {t2})."
                    f" Valid types are: {self.types.keys()}."
                )
            sigma, epsilon = self.types[(t1, t2)]
        c6 = 4 * epsilon * (sigma ** 6)
        c12 = 4 * epsilon * (sigma ** 12)
        self._force_obj.addBond(i, j, [qprod, c6, c12])

    _filters = {"pair"}

    def _additional_save(self, f):
        f.dump(self.types)

    def _additional_load(self, f):
        self.types = f.load()
