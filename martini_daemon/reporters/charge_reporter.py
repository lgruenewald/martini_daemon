from .reporter import Reporter
from ..utils import backup_try
import numpy as np


class ChargeReporter(Reporter):
    """
        Reporter that reports on the charge
    """

    def __init__(self):
        pass

    def on_set_xtc_path(self, xtc_name):
        filename = xtc_name + ".charges"
        backup_try(filename)

    def on_xtc_frame(self, pos, box, xtc_name):
        filename = xtc_name + ".charges"
        n = self._sysstar.len_particles()

        charges = []
        for i in range(n):
            _, charge, _ = self._sysstar.get_particle_details(i)
            charges.append(charge)
        charges_nd = np.array(charges, np.int32)
        with open(filename, "ab") as file:
            charges_nd.tofile(file)

