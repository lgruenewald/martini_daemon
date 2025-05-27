from .reporter import Reporter, read_compressed
import numpy as np
import struct

import zlib
from ..utils import backup_try, cross_box


def _append(arr, index, val, max, empty_val):
    """
        helper for simulated variable length numpy arrays
    """
    for i in range(max):
        if arr[index][i] == empty_val:
            arr[index][i] = val
            return
    assert False, f"More than {max} bonds for atom {index+1}."


def _has(arr, index, val, max, empty_val):
    """
        helper for simulated variable length numpy arrays
    """
    for i in range(max):
        if arr[index][i] == val:
            return True
        elif arr[index][i] == empty_val:
            return False
    return False


class BondReporter(Reporter):
    """
    Bond, constraint and vsite network graph connectivity reporter.
    In case of vsites, the virtual particle is considered bonded to all
    constructing particles.

    Can be used for:
    - analysis of reactions during a trajectory
    - visualization (if excluding bonds that cross pbc)
    - pbc whole or similar graph based trajectory manipulation
    Note: reports all bonds, including those crossing the PBC,
    reports at the same frames/frequency as the .xtc.

    constructor arguments:
    max_bonds_per_atom - maximum number of connections per atom, increase
    as needed, the default 12 should work for many use cases.
    max_atoms - only the first max_atoms atoms will be reported on (0=all)

    File format - zlib compressed binary data:
    - header:
        - n_atoms, 32 bit unsigned integer
        - max_bonds_per_atom, 32 bit unsigned integer
    - bond info - for each frame:
        - n_atoms "bonded to" connectivity data, array of max_bonds_per_atom
        32 bit unsigned integers
    Note: bonds are bidirectional - a bond between i and j means there is an
    entry at both i and j to the other.
    Note: assumes at most 2^32 - 1 atoms. Endianness used is native.
    """

    def __init__(self, max_bonds_per_atom=12, max_atoms=0):
        self.max_bonds_per_atom = max_bonds_per_atom
        assert self.max_bonds_per_atom > 0
        self.max_atoms = max_atoms

    def on_set_xtc_path(self, xtc_name):
        self._open_compressed(xtc_name + ".bonds")
        n = self._sysstar.len_particles()
        assert n > 0
        assert self.max_atoms < n, f"max_atoms ({self.max_atoms}) is larger than the total number of atoms ({n})."
        if self.max_atoms > 0:
            n = self.max_atoms
        self._write(struct.pack("=L", n))
        self._write(struct.pack("=L", self.max_bonds_per_atom))
        self.n = n

    def on_xtc_frame(self, i, pos, box, xtc_name):
        bonds = np.zeros((self.n, self.max_bonds_per_atom), dtype=np.uint32) - 1
        empty_val = 0xffffffff
        for force in self._sysstar.modular_forces:
            if force.is_instance("bond"):
                for id in range(len(force)):
                    i, j = force.get_members(id)
                    if i == j or i >= self.n or j >= self.n:
                        continue
                    if _has(bonds, i, j, self.max_bonds_per_atom, empty_val):
                        # multiple bonds between i and j? only add once
                        continue
                    _append(bonds, i, j, self.max_bonds_per_atom, empty_val)
                    _append(bonds, j, i, self.max_bonds_per_atom, empty_val)
            elif force.is_instance("vsite"):
                for id in range(len(force)):
                    vid, *others = force.get_members(id)
                    for other in others:
                        if vid == other or vid >= self.n or other >= self.n:
                            continue
                        if _has(bonds, vid, other, self.max_bonds_per_atom, empty_val):
                            continue
                        # multiple bonds between i and j? skip
                        _append(
                            bonds, vid, other, self.max_bonds_per_atom, empty_val
                        )
                        _append(
                            bonds, other, vid, self.max_bonds_per_atom, empty_val
                        )

        self._write(bonds.tobytes())


def read_bonds(path: str) -> tuple[int, int, int, int, np.array]:
    """
        Reads a file written by BondReporter.
        Returns:
        n_frames, n_atoms, max_bonds, empty_value, array
        where the array is of dimension (n_frames, n_atoms, max_bonds) and
        empty_value is used in the array to indicate non values.
    """
    data = read_compressed(path)
    n, max_bonds = struct.unpack("=LL", data[0:8])
    n_frames = (len(data) - 8) // (max_bonds * n * 4)
    assert (len(data) - 8) % (max_bonds * n * 4) == 0, "Wrong file format or incomplete file."
    array = np.frombuffer(data[8:], dtype=np.uint32)
    return (
        n_frames, n, max_bonds, 0xffffffff,
        array.reshape((n_frames, n, max_bonds))
    )


class VMDBondReporter(Reporter):
    """
    Reporter that prints a live list of bonds to a file, that can be
    visualized in VMD using daemon.tcl

    Limitation right now:
    Does not report bonds that cross the periodic boundary conditions,
    once a daemon aware pbc whole is written perhaps this can change"""

    handle = None

    def __init__(self):
        pass

    def on_set_xtc_path(self, xtc_name):
        self.bonds_name = xtc_name + ".bonds"
        backup_try(self.bonds_name)
        self.comp_obj = zlib.compressobj(6)
        self.handle = open(self.bonds_name, "wb")

    def on_xtc_frame(self, i, pos, box, _):
        bonds = "{" + "} {".join(
            map(
                lambda part:  # part[0] - i, part[1] - list of interactions
                " ".join(
                    map(
                        lambda ij:  # ij[0] - i, ij[1] - j
                            str(ij[1]) if ij[0] == part[0] else str(ij[0]),
                        map(
                            lambda inter: inter.get_members(),
                            filter(
                                lambda inter:
                                    not cross_box(
                                        pos[inter.get_members()[0]],
                                        pos[inter.get_members()[1]],
                                        box
                                    ),
                                filter(
                                    lambda inter: inter.is_instance("bond"),
                                    part[1]
                                )
                            )
                        )
                    )
                ),
                enumerate(self._topstar.interaction_list)
            )
        ) + "} "

        flat = self.comp_obj.compress(bonds.encode("utf-8"))

        self.handle.write(flat)

    def __del__(self):
        if self.handle is not None:
            self.handle.write(self.comp_obj.flush())
            self.handle.close()
