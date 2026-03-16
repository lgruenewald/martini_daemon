from .reporter import Reporter, read_compressed
import numpy as np
import struct


class AtomReporter(Reporter):
    """
        Reporter that reports the name, LJ type, mass and charge of all
        atoms in all xtc frames.

        file format - zlib compressed following:
        header: n_atoms (64 bit unsigned int)
        body - for each frame, without any padding or separator:
            - n_atoms s8 atom names
            - n_atoms s8 atom types
            - n_atoms f32 charges
            - n_atoms f32 masses
        f32 being a 32 bit float, s8 being an utf-8 string cut off at 8 bytes,
        null terminated if shorter than 8 bytes.

        note: endianness is the native endianness
    """

    def __init__(self):
        pass

    def on_set_xtc_path(self, xtc_name):
        self._open_compressed(xtc_name + ".atoms")
        n = self._sysstar.len_atoms()
        self._write(struct.pack("=Q", n))

    def on_xtc_frame(self, i, pos, box, xtc_name):
        n = self._sysstar.len_atoms()
        names = np.empty(n, dtype="a8")
        types = np.empty(n, dtype="a8")
        charges = np.empty(n, dtype=np.float32)
        masses = np.empty(n, dtype=np.float32)
        for i in range(n):
            name = self._sysstar.get_atom_name(i)
            type, charge, mass = self._sysstar.get_atom_details(i)
            # this will silently cut off everything after 8 bytes
            names[i] = name.encode("utf-8")
            types[i] = type.encode("utf-8")
            charges[i] = charge
            masses[i] = mass
        self._write(names.tobytes())
        self._write(types.tobytes())
        self._write(charges.tobytes())
        self._write(masses.tobytes())


def read_atoms(path):
    """
        Reads a file created by this reporter.

        Returns:
        n_frames, n_atoms, names, types, charges, masses
        where the last 4 are numpy arrays of dimension (n_frames, n_atoms)
    """
    data = read_compressed(path)
    n_atoms, = struct.unpack("=Q", data[0:8])
    # length for a single atom: 8+8+4+4 = 24 bytes
    n_frames = (len(data) - 8) // (n_atoms * 24)
    assert (len(data) - 8) % (n_atoms * 24) == 0, "The .atoms file is incomplete or the wrong file type."
    names = np.empty((n_frames, n_atoms), dtype="a8")
    types = np.empty((n_frames, n_atoms), dtype="a8")
    charges = np.empty((n_frames, n_atoms), dtype=np.float32)
    masses = np.empty((n_frames, n_atoms), dtype=np.float32)
    for frame in range(n_frames):
        start = 8 + (frame * 24 * n_atoms)
        names[frame] = np.frombuffer(
            data[start:start+(8*n_atoms)], dtype="a8"
        )
        types[frame] = np.frombuffer(
            data[start+8*n_atoms:start+16*n_atoms], dtype="a8"
        )
        charges[frame] = np.frombuffer(
            data[start+16*n_atoms:start+20*n_atoms], dtype=np.float32
        )
        masses[frame] = np.frombuffer(
            data[start+20*n_atoms:start+24*n_atoms], dtype=np.float32
        )
    return n_frames, n_atoms, names, types, charges, masses
