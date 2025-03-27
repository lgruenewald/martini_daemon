from .reporter import Reporter
from ..utils import backup_try
import numpy as np


class AllBondReporter(Reporter):
    """Reporter that prints a live list of bonds to a file, that can be
    visualized in VMD using daemon.tcl, including pbc crossing"""

    def on_set_xtc_path(self, xtc_name):
        bonds_name = xtc_name + "_all.bonds"
        backup_try(bonds_name)

    def on_xtc_frame(self, pos, box, xtc_name):
        bonds_name = xtc_name + "_all.bonds"
        n = self._sysstar.len_particles()

        bonds = [[] for _ in range(n)]
        for force in self._sysstar.modular_forces:
            if not force.is_instance("bond"):
                continue
            for entry in filter(None, force._list):
                i, j, *_ = entry
                bonds[i].append(j)
                bonds[j].append(i)
        # put a -1 between every list, which is how daemon.tcl reads it
        for c in range(len(bonds)):
            bonds[c].append(-1)
        flat = [x for xs in bonds for x in xs]
        flat = np.array(flat, np.int32)
        with open(bonds_name, "ab") as file:
            flat.tofile(file)

