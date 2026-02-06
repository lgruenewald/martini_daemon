from .reporter import Reporter, read_compressed
import numpy as np
import struct


class SysStarDump(Reporter):
    """
        Dumps all info from S* -- all atom details and all interactions
        in a human readable plaintext file. Dumps it at the start and
        when there is any change.
    """

    def __init__(self):
        pass

    def on_set_xtc_path(self, xtc_name):
        self._open(xtc_name + ".sstar")
        n = self._sysstar.len_atoms()
        assert n > 0
        self._print("SysStar dump")
        self._print(f"# of atoms {n}")
        self._print("")
        self._write_frame(0)

    def post_modification(self, i, name, reactions):
        self._write_frame(i)

    def _write_frame(self, i):
        self._print(f"==== Frame {i} ====")

        sys = self._sysstar
        self._print("Atom list")
        self._print("(name, resname, type, charge, mass, sc_lam, sc_alpha)")
        for deets in sys._atom_list:
            self._print(f"{deets}")

        self._print("")
        self._print("Force list")
        for force in sys.modular_forces:
            self._print(f"Force {force}")
            if force._list is None:
                self._print("List is none")
            else:
                for deets in force._list:
                    self._print(f"{deets}")


        self._print(f"==== End of Frame {i} ====")
        self._print("")
