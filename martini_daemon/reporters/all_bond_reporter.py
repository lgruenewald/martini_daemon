from .reporter import Reporter
from ..utils import backup_try
import numpy as np
import zlib


class AllBondReporter(Reporter):
    """Reporter for python analysis scripts"""

    handle = None
    def __init__(self):
        pass

    def on_set_xtc_path(self, xtc_name):
        self.bonds_name = xtc_name + ".all_bonds"
        backup_try(self.bonds_name)
        self.comp_obj = zlib.compressobj(6)
        self.handle = open(self.bonds_name, "wb")

    def on_xtc_frame(self, pos, box, xtc_name):
        n = self._sysstar.len_particles()

        bonds = np.zeros((n, 12), dtype=np.int32) - 1
        for part_id, part in enumerate(self._topstar.interaction_list):
            bi = 0
            for inter in part:
                if not inter.is_instance("bond"):
                    continue
                i, j = inter.get_members()
                bonds[part_id][bi] = i if part_id == j else j
                bi += 1
                if bi > 12:
                    raise ValueError("More than 12 bonds per atom")

        self.handle.write(self.comp_obj.compress(bonds.tobytes()))

    def __del__(self):
        if self.handle is not None:
            self.handle.write(self.comp_obj.flush())
            self.handle.close()
