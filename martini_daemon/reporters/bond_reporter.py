from .reporter import Reporter, read_compressed
import numpy as np
import struct

def collect_bonds(n, sstar, constraint_only=False):
    bonds_len = 0
    vsite_len = 0
    bond_filter = "constraint" if constraint_only else "bond"
    for force in sstar.modular_forces:
        if force.is_instance(bond_filter):
            for id in range(len(force)):
                if force.get_members(id) is not None:
                    bonds_len += 1
        elif force.is_instance("vsite"):
            for id in range(len(force)):
                if force.get_members(id) is not None:
                    vsite_len += len(force.get_members(id)) - 1
    # allocate
    n_bonds = bonds_len + vsite_len
    bonds = np.empty((n_bonds, 2), dtype=np.uint32)
    # write
    bond_index = 0
    for force in sstar.modular_forces:
        if force.is_instance(bond_filter):
            for id in range(len(force)):
                members = force.get_members(id)
                if members is None:
                    continue
                i, j = members
                if i == j or i >= n or j >= n:
                    continue
                bonds[bond_index][0] = i
                bonds[bond_index][1] = j
                bond_index += 1
        elif force.is_instance("vsite"):
            for id in range(len(force)):
                members = force.get_members(id)
                if members is None:
                    continue
                vid, *others = members
                for other in others:
                    if vid == other or vid >= n or other >= n:
                        continue
                    bonds[bond_index][0] = vid
                    bonds[bond_index][1] = other
                    bond_index += 1

    assert bond_index == n_bonds
    return n_bonds, bonds

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

    def __init__(self, max_atoms=0):
        self.max_atoms = max_atoms

    def on_set_xtc_path(self, xtc_name):
        self._open_compressed(xtc_name + ".bonds")
        n = self._sysstar.len_atoms()
        assert n > 0
        assert self.max_atoms < n, (
            f"max_atoms ({self.max_atoms}) is larger "
            f"than the total number of atoms ({n})."
        )
        if self.max_atoms > 0:
            n = self.max_atoms
        self.n = n
        self._write(struct.pack("=Q", n))

    def on_xtc_frame(self, frame_index, pos, box, xtc_name):
        n_bonds, bonds = collect_bonds(self.n, self._sysstar)
        self._write(struct.pack("=Q", n_bonds))
        self._write(bonds.tobytes())


def read_bonds(path: str, read_n_atoms=True) -> list[np.ndarray]:
    """
        Reads a file written by BondReporter.
        Returns: n_frames, n_atoms, frames
        frames = a list of frames, each frame containing a numpy
        array of (n_bonds, 2) shape, where n_bonds can vary per frame.

        set read_n_atoms to False when reading old bond reporter outputs.
    """
    data = read_compressed(path)
    if read_n_atoms:
        n_atoms, = struct.unpack("=Q", data[0:8])
        frame_start = 8
    else:
        n_atoms = 0
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
    n_frames = len(frames)
    return n_frames, n_atoms, frames
