from .force import Force
import openmm as mm  # type: ignore[import-untyped]
from openmm.unit import nanometer  # type: ignore[import-untyped]


class Pairs(Force):
    """Pairs interactions in the system
    values are (i: part_id, j: part_id, type, qprod, param1, param2)

    if type is 0 use pair types, param1 and param2 are LJ types

    if type is 1 param1 and param2 are c6 and c12"""

    types: dict[tuple[str, str], tuple[float, float]]

    def __init__(self, sysstar):
        super().__init__(sysstar)
        self.types = {}

    def add_type(self, type1, type2, sigma, epsilon):
        # TODO move V,W conversion to top_parser and unify in single helper
        c6 = 4 * epsilon * (sigma ** 6)
        c12 = 4 * epsilon * (sigma ** 12)
        self.types[(type1, type2)] = (c6, c12)
        self.types[(type2, type1)] = (c6, c12)

    def _addbond(self, i, j, p1, p2):
        # p1, p2 are c6 and c12 in this order if type is 1
        # p1, p2 are atom types if type is 0
        t1, q1, _ = self._sysstar.get_particle_details(i)
        t2, q2, _ = self._sysstar.get_particle_details(j)
        qprod = q1 * q2
        if p1 is None or p2 is None:
            p1, p2 = self.types[(t1, t2)]
        self._force_obj.addBond(i, j, [qprod, p1, p2])

    def _build(self):
        epsilon_r = self._sysstar.epsilon_r
        self._force_obj = mm.CustomBondForce(
            "LJ + ES;"
            "LJ = (C12 / r^12 - C6 / r^6);"
            "corr = (C12 / rcut^12 - C6 / rcut^6);"
            "ES = f*qprod/epsilon_r/r;"
            f"epsilon_r = {epsilon_r};"
            "f = 138.935458;"
        )
        self._force_obj.addPerBondParameter("qprod")
        self._force_obj.addPerBondParameter("C6")
        self._force_obj.addPerBondParameter("C12")

        for (i, j, p1, p2) in filter(None, self._list):
            self._addbond(i, j, p1, p2)

    def add(self, members, params):
        # process params into something storeable
        i, j = members
        # default behavior is to use the [pairtypes] directive's values
        p1, p2 = None, None
        if len(params) == 2:
            sigma, epsilon = params
            p1 = 4 * epsilon * (sigma ** 6)  # C6
            p2 = 4 * epsilon * (sigma ** 12)  # C12

        # typical add() boilerplate
        self._list.append((i, j, p1, p2))
        if not self._rebuild:
            self._addbond(i, j, p1, p2)
            self._sysstar._reinitialize = True
        return self._interaction()

    def get_members(self, id):
        i, j, _, _ = self._list[id]
        return [i, j]

    def update_params(self, id, length, kb, kcub):
        raise NotImplementedError

    def is_instance(self, filter):
        return filter == "pairs"

