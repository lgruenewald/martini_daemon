from .reporter import Reporter
from ..utils import backup_try
import numpy as np


class TypeReporter(Reporter):
    """
        Reporter that reports on the LJ types (with numerical indices)
    """

    def pre_steps(self, xtc_name):
        bonds_name = xtc_name + "_types.npy"
        backup_try(bonds_name)

    def step(self, pos, box, xtc_name):
        bonds_name = xtc_name + "_types.npy"
        n = len(self._sysstar._part_list)

        types = []
        for i in range(n):
            _, type, _, _ = self._sysstar.get_particle_details(i)
            type_id = self._sysstar.nonbonded_force.use_atom_type(type)
            types.append(type_id)
        types = np.array(types, np.int32)
        with open(bonds_name, "ab") as file:
            types.tofile(file)

