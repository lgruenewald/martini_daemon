from .reporter import Reporter
from ..utils import backup_try
import numpy as np


class TypeReporter(Reporter):
    """
        Reporter that reports on the LJ types (with numerical indices)
    """

    def __init__(self):
        pass

    def on_set_xtc_path(self, xtc_name):
        bonds_name = xtc_name + "_types.npy"
        backup_try(bonds_name)

    def on_xtc_frame(self, i, pos, box, xtc_name):
        bonds_name = xtc_name + "_types.npy"
        n = self._sysstar.len_particles()

        types = []
        for i in range(n):
            type, _, _ = self._sysstar.get_particle_details(i)
            type_id = self._sysstar.nonbonded_force.use_atom_type(type)
            types.append(type_id)
        types_nd = np.array(types, np.int32)
        with open(bonds_name, "ab") as file:
            types_nd.tofile(file)

