from .force import Force
import openmm as mm


class Cmap(Force):
    _members = 5

    def __init__(self, sysstar):
        super().__init__(sysstar)
        self.types = {}
        self.maps = []

    def add_type(self, members, params):
        # rearrangement as in
        # https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/gromacstopfile.py
        # lines 959-984
        size, p2, *ps = params
        assert size == p2
        assert len(ps) == size*size
        cmap = []
        midpoint = size // 2

        # based on testing:
        # gromacs CMAPs start at -180,-180
        # rows are the second dihedral, going from -180 to +180
        # columns are the first dihedral, going from -180 to +180

        for n in range(size):
            column = (n + midpoint) % size
            for o in range(size):
                row = (o + midpoint) % size
                cmap.append(ps[size * row + column])

        self.maps.append((size, cmap))
        if self._force_obj is not None:
            i = self._force_obj.addMap(size, cmap),
            assert len(self.maps) - 1 == i
        self.types[tuple(members)] = len(self.maps) - 1

    def _set_force_obj(self):
        self._force_obj = mm.CMAPTorsionForce()
        for index, (size, cmap) in enumerate(self.maps):
            i = self._force_obj.addMap(size, cmap)
            assert i == index

    def _add_to_force_obj(self, params):
        i, j, k, l, m = params
        types = tuple([
            self._sysstar.get_atom_details(x)[0]
            for x in [i, j, k, l, m]
        ])
        default = self.types.get(types)
        if default is None:
            raise ValueError(f"Unknown CMAP type for {types}.")

        self._force_obj.addTorsion(
            default,
            i, j, k, l,
            j, k, l, m
        )

    _filters = {"cmap"}

    def _additional_save(self, f):
        f.dump(self.types)
        f.dump(self.maps)

    def _additional_load(self, f):
        self.types = f.load()
        self.maps = f.load()
