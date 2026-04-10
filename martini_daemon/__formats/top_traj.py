import struct
import zlib
from collections.abc import Collection
from dataclasses import dataclass

from ..__rust import BondGraph

# TODO improve docs
# TODO initial molecules should contain number of atoms per molecule as well
# TODO reaction information
# TODO fragment information per frame?


class TopTrajWriter:
    """The Topology-Trajectory File Format .toptraj is described here.

    Goals of this format:
    - Store per-trajectory data for:
        - Simulation name
        - LIST OF INITIAL MOLECULES - a list of molecule names and counts. Only one per trajectory. Stored because
          it can be of interest for some analysis scripts.
    - Store per-frame data for:
        - Frame number, simulation step, simulation time, number of atoms
        - PER-ATOM DATA: atom names, residue names, residue IDs, atom types, charges, masses
        - LIST OF BONDS - all bonds, constraints, virtual sites. No data on the bond strengths, as the connectivity
          network is the most important, bond strength is a force field detail. Virtual sites always have bonds between
          the constructing atoms and the virtual particle, no extra bonds added among constructing atoms.
          It also sets a maximum of one bond between each pair of atoms, as the number of bonds as defined in a force
          field topology may not be so relevant for a connectivity graph.
    - Optimize disk space and sequential access (both reading and writing)
    - No forward compatibility, not even across minor version
        - it is a non-goal for old reader software to be able to read newer files
        - readers should still be able to read older files thanks to the version string

    Assumptions made about the simulation:
        - less than 2^32 frames in total (the MD step can exceed 2^32 though)
        - less than 2^32 atoms
        - less than 2^32 residues

    Description of the binary layout:

    - TopTraj files consist of a header and body.
    - The header contains metadata only, and is present as raw bytes.
    - The body contains all the info described above. The header is compressed using ZLIB. ZLIB writes its own tiny
    header at the start of the body which contains the compression details. A compression level of 6 is recommended.
    - integers are all little endian and unsigned.
    - all strings are UTF-8 encoded and encoded as pascal strings
        - One byte length, followed by the bytes of the string.
        - Writers may truncate as they wish, even resulting in
        incomplete UTF-8 codepoints.
        - Readers should stop parsing UTF-8 at the invalid
        UTF-8 character, and skip the remainder of the bytes.
    - floats are 32 bit IEEE floats.
    - doubles are 64 bit IEEE floats.

    The header contains the following information:
    - The 6 byte Magic number "0xc0TOPTR"
    - 1 byte for major version of format (currently 1)
    - 1 byte for minor version of format (currently 0)

    The body contains the following information:
    - Simulation name
        - string - simulation title
    - List of initial molecules:
        - 4 byte integer - number of entries
        - List of entries:
            - string - molecule name
            - 4 byte integer - number of molecules in this entry
    - CRC32 checksum of header bytes+simulation name+list of initial molecules (as raw uncompressed bytes for the latter two)
    - List of frames:
        - 4 byte integer - frame number, must be one greater than previous frame, must be 0 for the first frame
        - 4 byte integer - number of atoms (n_atoms)
        - 8 byte integer - simulation step
        - double - simulation time, in picoseconds
        - List of atom names
            - n_atoms strings
        - Similar setup for resname, resid, atom type, charges, mass. Payload type:
            - resname - string
            - resid - 4 byte integer
            - atom type - string
            - charge - float
            - mass - float
        - List of bonds:
            - in no particular order, however guaranteed to have at most one bond per i-j combination and no i-i
            - 8 byte integer n_bonds
            - 4 byte integer, index of atom i
            - 4 byte integer, index of atom j
        - CRC32 checksum of frame data
    """

    def __init__(self, path: str, title: str, initial_molecules: list[tuple[str, int]]):
        """Creates a Topology Trajectory writer.

        Note: this class owns a file handle. If using it directly, call .finish() manually when done!
        """
        # version 1.0
        header = b"\xc0TOPTR\x01\x00"
        title_bytes = title.encode("utf-8")
        title_len = min(255, len(title_bytes))
        n_initial = len(initial_molecules)
        body = struct.pack(f"<B{title_len}sI", title_len, title_bytes, n_initial)
        for name, count in initial_molecules:
            name_bytes = name.encode("utf-8")
            name_len = min(255, len(name_bytes))
            body += struct.pack(
                # TODO test with 255< bytes
                f"<B{name_len}sI",
                name_len,
                name_bytes,
                count,
            )
        checksum = zlib.crc32(header + body)
        body += struct.pack("<I", checksum)

        self.handle = open(path, "wb")  # noqa: SIM115
        self.handle.write(header)
        self.comp = zlib.compressobj()
        self.handle.write(self.comp.compress(body))
        self.handle.write(self.comp.flush(zlib.Z_PARTIAL_FLUSH))
        self.handle.flush()
        self.last_frame = -1
        self.crc32 = 0
        # for now, this keeps a lot of information in memory both here
        # and in system
        #
        # if this extra memory consumption becomes a problem this should be
        # more closely integrated within System
        #
        # currently it's put here during prototyping to keep it modular
        self.frame = None
        self.previous_frame = None

    def new_frame(self, frame_num: int, sim_step: int, time_ps: float, n_atoms: int):
        """Create a new frame."""
        assert frame_num - self.last_frame == 1, (
            "Frames passed to TopTrajWriter must be in a sequence. "
            f"Got frame num {frame_num}, expected {self.last_frame + 1}."
        )
        assert self.frame is None, "Only call new_frame after write_frame()!"
        self.last_frame = frame_num

        self.write(struct.pack("<IIQd", frame_num, n_atoms, sim_step, time_ps))

        self.frame = {"n_atoms": n_atoms}

    def register_frame_atoms(
        self,
        names: Collection[str],
        res_names: Collection[str],
        res_ids: Collection[int],
        atom_types: Collection[str],
        charges: Collection[float],
        masses: Collection[float],
    ):
        """Writes the current frame atom information to disk."""
        if self.frame is None:
            raise ValueError("Must call new_frame() first!")

        if self.frame.get("atoms") is not None:
            raise ValueError("Should only call register_frame_atoms once per frame!")

        assert (
            len(names)
            == len(res_names)
            == len(res_ids)
            == len(atom_types)
            == len(charges)
            == len(masses)
        )
        assert len(names) == self.frame["n_atoms"]

        for name in names:
            name_bytes = name.encode("utf-8")
            self.write(
                struct.pack(f"<B{len(name_bytes)}s", len(name_bytes), name_bytes)
            )

        for name in res_names:
            name_bytes = name.encode("utf-8")
            self.write(
                struct.pack(f"<B{len(name_bytes)}s", len(name_bytes), name_bytes)
            )

        for res_id in res_ids:
            self.write(struct.pack("<I", res_id))

        for atom_type in atom_types:
            atom_type_bytes = atom_type.encode("utf-8")
            self.write(
                struct.pack(
                    f"<B{len(atom_type_bytes)}s", len(atom_type_bytes), atom_type_bytes
                )
            )

        for charge in charges:
            self.write(struct.pack("<f", charge))

        for mass in masses:
            self.write(struct.pack("<f", mass))

        self.frame["atoms"] = True

    def register_frame_bonds(self, bonds: BondGraph) -> None:
        """Writes the bonds for the current frame to disk."""
        if self.frame is None:
            raise ValueError("Must call new_frame() first!")
        if self.frame.get("bonds") is not None:
            raise ValueError("Must only call register_frame_bonds once per frame!")

        self.frame["bonds"] = True
        raw_bytes = bytearray()

        bonds_len = len(raw_bytes)
        raw_bytes += struct.pack("<Q", 0)
        n_bonds = 0
        wrote = set()
        for i, j in bonds.to_list():
            if i == j:
                continue

            smaller = min(i, j)
            larger = max(i, j)
            if (smaller, larger) in wrote:
                continue
            wrote.add((smaller, larger))
            n_bonds += 1
            raw_bytes += struct.pack("<II", smaller, larger)
        raw_bytes[bonds_len : bonds_len + 8] = struct.pack("<Q", n_bonds)
        self.write(raw_bytes)

    def write_frame(self):
        """Finishes writing the current frame to disk. Must call register_frame_atoms and register_frame_bonds exactly once first."""
        self.frame = None
        self.write_crc32()
        self.flush(partial=True)

    def write(self, raw_bytes: bytes | bytearray) -> None:
        assert self.handle is not None
        assert self.comp is not None
        self.handle.write(self.comp.compress(raw_bytes))
        self.crc32 = zlib.crc32(raw_bytes, self.crc32)

    def write_crc32(self):
        self.write(struct.pack("<I", self.crc32))
        self.crc32 = 0

    def flush(self, partial=False):
        assert self.handle is not None
        assert self.comp is not None
        if partial:
            self.handle.write(self.comp.flush(zlib.Z_PARTIAL_FLUSH))
        else:
            self.handle.write(self.comp.flush())
        self.handle.flush()

    def finish(self):
        assert self.handle is not None
        self.flush()
        self.handle.close()
        self.handle = None
        self.comp = None


@dataclass
class TopTrajFrame:
    frame_index: int
    sim_step: int
    sim_time: float
    n_atoms: int
    names: list[str]
    res_names: list[str]
    res_ids: list[int]
    atom_types: list[str]
    charges: list[float]
    masses: list[float]
    bonds: list[tuple[int, int]]


class TopTrajReader:
    """For a description of the file format, see TopTrajWriter."""

    def __init__(self, path: str):
        self.path = path
        with open(path, "rb") as f:
            self.header = f.read(8)
            assert self.header == b"\xc0TOPTR\x01\x00"
            self.content = zlib.decompress(f.read())
        title_len = self.content[0]
        (self.title,) = struct.unpack(f"<{title_len}s", self.content[1 : 1 + title_len])
        self.i = 1 + title_len
        (n_init,) = struct.unpack("<I", self.content[self.i : self.i + 4])
        self.i += 4
        self.initial_molecules = []
        for i in range(n_init):
            name_len = self.content[self.i]
            self.i += 1
            (name,) = struct.unpack(
                f"<{name_len}s", self.content[self.i : self.i + name_len]
            )
            self.i += name_len
            (n,) = struct.unpack("<I", self.content[self.i : self.i + 4])
            self.i += 4
            self.initial_molecules.append((name, n))

        chunk1 = self.header + self.content[: self.i]
        (crc,) = struct.unpack("<I", self.content[self.i : self.i + 4])
        assert zlib.crc32(chunk1) == crc, (
            f"File {path} appears to be corrupt. Header CRC32 {crc} doesn't match {zlib.crc32(chunk1)}."
        )
        self.i += 4
        self.frame = 0

    def __read_strings(self, n):
        res = []
        for i in range(n):
            len_ = self.content[self.i]
            self.i += 1
            res.append(
                struct.unpack(f"<{len_}s", self.content[self.i : self.i + len_])[0]
            )
            self.i += len_
        return res

    def __read_any(self, n, byte_format, n_bytes):
        res = []
        for i in range(n):
            res.append(
                struct.unpack(byte_format, self.content[self.i : self.i + n_bytes])[0]
            )
            self.i += n_bytes
        return res

    def read_frame(self) -> TopTrajFrame | None:
        if self.i >= len(self.content):
            return None

        frame_start = self.i
        frame, n_atoms, sim_step, sim_time = struct.unpack(
            "<IIQd", self.content[self.i : self.i + 24]
        )
        self.i += 24
        names = self.__read_strings(n_atoms)
        res_names = self.__read_strings(n_atoms)
        res_ids = self.__read_any(n_atoms, "<I", 4)
        types = self.__read_strings(n_atoms)
        charges = self.__read_any(n_atoms, "<f", 4)
        masses = self.__read_any(n_atoms, "<f", 4)
        (n_bonds,) = struct.unpack("<Q", self.content[self.i : self.i + 8])
        bonds = []
        self.i += 8
        for i in range(n_bonds):
            bonds.append(struct.unpack("<II", self.content[self.i : self.i + 8]))
            self.i += 8
        (crc,) = struct.unpack("<I", self.content[self.i : self.i + 4])
        self.i += 4
        assert zlib.crc32(self.content[frame_start : self.i - 4]) == crc, (
            f"Frame {frame} corrupt, CRC32 mismatch."
        )

        return TopTrajFrame(
            frame,
            sim_step,
            sim_time,
            n_atoms,
            names,
            res_names,
            res_ids,
            types,
            charges,
            masses,
            bonds,
        )
