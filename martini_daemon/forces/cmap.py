from .force import Force
import openmm as mm


class Cmap(Force):
    _members = 5

    def __init__(self, sysstar):
        super().__init__(sysstar)
        self.types = {}

    def add_type(self, members, params):
        self.types[tuple(members)] = tuple(params)

    def _set_force_obj(self):
        self._force_obj = mm.CMAPTorsionForce()

    def _add_to_force_obj(self, params):
        i, j, k, L, m = params
        types = tuple([
            self._sysstar.get_atom_details(x)[0]
            for x in [i, j, k, L, m]
        ])
        default = self.types.get(types)
        if default is None:
            raise ValueError(f"Unknown CMAP type for {types}.")
        p1, p2, *ps = default
        assert p1 == p2
        
        # rearrangement as in https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/gromacstopfile.py
        # lines 959-984
        map = []
        size = p1
        midpoint = size // 2
        for n in range(size):
            column = (n + midpoint) % size
            for o in range(size):
                row = (o + midpoint) % size
                map.append(ps[size * row + column])

        self._force_obj.addTorsion(
            self._force_obj.addMap(size, map),
            i, j, k, L,
            j, k, L, m
        )

    _filters = {"cmap"}

    def _additional_save(self, f):
        f.dump(self.types)

    def _additional_load(self, f):
        self.types = f.load()
