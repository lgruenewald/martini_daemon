from .reporter import Reporter
from ..utils import backup_try, cross_box
import numpy as np


class BondReporter(Reporter):
    """Reporter that prints a live list of bonds to a file, that can be
    visualized in VMD using daemon.tcl

    Limitation right now:
    Does not report bonds that cross the periodic boundary conditions,
    once a daemon aware pbc whole is written perhaps this can change"""

    def pre_steps(self, xtc_name):
        bonds_name = xtc_name + "_bonds.npy"
        backup_try(bonds_name)

    def step(self, pos, box, xtc_name):
        bonds_name = xtc_name + "_bonds.npy"
        n = self._sysstar.len_particles()

        bonds = [[] for _ in range(n)]
        for force in self._sysstar.modular_forces:
            if not force.visualize_as_bond:
                continue
            for entry in filter(None, force._list):
                i, j, *_ = entry
                if cross_box(pos[i], pos[j], box):
                    continue
                bonds[i].append(j)
                bonds[j].append(i)
        # put a -1 between every list, which is how daemon.tcl reads it
        for c in range(len(bonds)):
            bonds[c].append(-1)
        flat = [x for xs in bonds for x in xs]
        flat = np.array(flat, np.int32)
        with open(bonds_name, "ab") as file:
            flat.tofile(file)

