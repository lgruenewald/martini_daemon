from .reporter import Reporter, read_compressed
import numpy as np
import struct


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
    reports at the same frames/frequency as the .xtc. Given multiple bonds
    between the same pair of particles, reports them multiple times.

    constructor arguments:
    max_atoms - only the first max_atoms atoms will be reported on (0=all)

    File format - zlib compressed binary data:
    - bond info - for each frame:
        - n_bonds - number of bonds in this frame (64 bit unsigned integer)
        - i, j - atom indices for a bond (both 32 bit unsigned integers)
    Note: assumes at most 2^32 atoms. Endianness used is native.
    """

    def __init__(self, max_bonds_per_atom=12, max_atoms=0):
        self.max_bonds_per_atom = max_bonds_per_atom
        assert self.max_bonds_per_atom > 0
        self.max_atoms = max_atoms

    def on_set_xtc_path(self, xtc_name):
        self._open_compressed(xtc_name + ".bonds")
        n = self._sysstar.len_atoms()
        assert n > 0
        assert self.max_atoms < n, f"max_atoms ({self.max_atoms}) is larger than the total number of atoms ({n})."
        if self.max_atoms > 0:
            n = self.max_atoms
        self.n = n

    def on_xtc_frame(self, frame_index, pos, box, xtc_name):
        bonds_len = sum(
            map(
                lambda force: len(force),
                filter(
                    lambda force: force.is_instance("bond"),
                    self._sysstar.modular_forces
                )
            )
        )
        vsite_len = 0
        for vsite in filter(lambda force: force.is_instance("vsite"), self._sysstar.modular_forces):
            for i in range(len(vsite)):
                # one vsite for every vid - constructing particle
                vsite_len += len(vsite.get_members(i)) - 1
        self._write(struct.pack("=Q", bonds_len + vsite_len))
        bonds = np.empty((bonds_len + vsite_len, 2), dtype=np.uint32)
        bond_index = 0
        for force in self._sysstar.modular_forces:
            if force.is_instance("bond"):
                for id in range(len(force)):
                    i, j = force.get_members(id)
                    if i == j or i >= self.n or j >= self.n:
                        continue
                    bonds[bond_index][0] = i
                    bonds[bond_index][1] = j
                    bond_index += 1
            elif force.is_instance("vsite"):
                for id in range(len(force)):
                    vid, *others = force.get_members(id)
                    for other in others:
                        if vid == other or vid >= self.n or other >= self.n:
                            continue
                        bonds[bond_index][0] = vid
                        bonds[bond_index][1] = other
                        bond_index += 1

        assert bond_index == vsite_len + bonds_len
        self._write(bonds.tobytes())


def read_bonds(path: str) -> tuple[int, int, int, int, np.array]:
    """
        Reads a file written by BondReporter.
        Returns: frames - a list of frames, each frame containing a numpy
        array of (n_bonds, 2) shape, where n_bonds can vary per frame.
    """
    data = read_compressed(path)
    frame_start = 0
    frames = []
    while frame_start < len(data):
        n_bonds, = struct.unpack("=Q", data[frame_start:frame_start+8])
        frame_end = frame_start + 8 + 8*n_bonds
        frames.append(
            np.frombuffer(data[frame_start+8:frame_end], dtype=np.uint32)
            .reshape((n_bonds, 2))
        )
        frame_start = frame_end
    return frames
