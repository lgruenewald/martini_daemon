from __future__ import annotations

import struct
import zlib
from collections.abc import Collection
from dataclasses import dataclass
from io import SEEK_CUR
from typing import Any

from ..__rust import BondGraph


class TopTrajWriter:
    """The Topology-Trajectory File Format .toptraj is described here.

    Goals of this format:
    - Store per-trajectory data in the header:
        - Simulation name
        - LIST OF INITIAL MOLECULES - a list of molecule names, counts and atoms per molecule.
          Only one per trajectory. Stored because it can be of interest for some analysis scripts.
        - residue names, residue IDs
    - Store per-frame data for:
        - Frame number, simulation step, simulation time, number of atoms
        - PER-ATOM DATA: atom names, atom types, charges, masses
        - LIST OF BONDS - all bonds, constraints, virtual sites. No data on the bond strengths, as the connectivity
          network is the most important, bond strength is a force field detail. Virtual sites always have bonds between
          the constructing atoms and the virtual particle, no extra bonds added among constructing atoms.
          It also sets a maximum of one bond between each pair of atoms, as the number of bonds as defined in a force
          field topology may not be so relevant for a connectivity graph.
    - Header and Per-frame compression using zlib, mmap-able file format for quick reading without loading it all
      into memory.
    - No forward compatibility, not even across minor version
        - it is a non-goal for old reader software to be able to read newer files
        - readers should still be able to read older files thanks to the version string

    Assumptions made about the simulation:
        - less than 2^32 frames in total (the MD step can exceed 2^32 though)
        - less than 2^32 atoms, a constant number of atoms per simulation
        - less than 2^32 residues

    Description of the binary layout:

    - TopTraj files consist of a magic number, header and body.
    - The header, and each frame is compressed using ZLIB. At the start, the size of the compressed chunk is given
      as an 8 byte integer. This is followed by the size of the (original) decompressed chunk, as an 8 byte integer.
      Then, ZLIB writes its own tiny header at the start of the body which contains the compression details.
      A compression level of 6 is recommended. At the end of each compressed chunk, a CRC32
      4 byte integer is written about the compressed contents. This crc32 is not compressed.
    - integers are all little endian and unsigned.
    - all strings are UTF-8 encoded and encoded as pascal strings
        - One byte length, followed by the bytes of the string.
        - Writers may truncate as they wish, even resulting in incomplete UTF-8 codepoints.
        - Readers should stop parsing UTF-8 at the invalid character, and skip the remainder of the bytes
        - There are no null terminating characters!
    - floats are 32 bit IEEE floats.
    - doubles are 64 bit IEEE floats.

    The magic number contains the following information:

    - The 6 byte Magic number "0xc0TOPTR"
    - 1 byte for major version of format (currently 1)
    - 1 byte for minor version of format (currently 0)

    The header contains the following information:

    - Simulation name
        - string - simulation title
    - List of initial molecules:
        - 4 byte integer - number of entries
        - List of entries:
            - string - molecule name
            - 4 byte integer - number of molecules in this entry
            - 4 byte integer - number of atoms per molecule
    - Initial n_atoms as a 4 byte integer
    - n_atoms strings for resnames
    - n_atoms 4 byte integers for resids

    Each frame contains the following information:

    - 4 byte integer - frame number, must be one greater than previous frame, must be 0 for the first frame
    - 4 byte integer - number of atoms (n_atoms)
    - 8 byte integer - simulation step
    - double - simulation time, in picoseconds
    - n_atoms strings for atom names
    - n_atoms strings for atom types
    - n_atoms floats for charge
    - n_atoms floats for mass
    - List of bonds:
        - in no particular order, however guaranteed to have at most one bond per i-j combination and no i-i
        - 8 byte integer n_bonds
        - 4 byte integer, index of atom i
        - 4 byte integer, index of atom j
    """

    def __init__(
        self,
        path: str,
        title: str,
        initial_molecules: list[tuple[str, int, int]],
        res_names: list[str],
        res_ids: list[int],
        append: bool = False,
        truncate: None | int = None,
    ) -> None:
        """Creates a Topology Trajectory writer.

        Note: this class owns a file handle. If using it directly, call .finish() manually when done!

        :param path: The path to write to.
        :param title: The simulation title to write to the file.
        :param initial_molecules: A list of tuples (name, count), corresponding to molecule names and counts at the
            start of the simulation.
        :param append: Whether to append to an existing file.
        :param truncate: If appending, the last MD step to keep.
        """
        self.__buffer = bytearray()
        self.__last_frame = -1
        self.__frame = None
        self.n_atoms = len(res_names)
        assert len(res_names) == len(res_ids)

        magic_1_0 = b"\xc0TOPTR\x01\x00"

        if not append:
            self.__handle = open(path, "wb")  # noqa: SIM115
            self.__handle.write(magic_1_0)

            title_bytes = title.encode("utf-8")
            title_len = min(255, len(title_bytes))
            n_initial = len(initial_molecules)
            self.__write(
                struct.pack(f"<B{title_len}sI", title_len, title_bytes, n_initial)
            )
            for name, count, atom_per_mol in initial_molecules:
                name_bytes = name.encode("utf-8")
                name_len = min(255, len(name_bytes))
                self.__write(
                    struct.pack(
                        f"<B{name_len}sII", name_len, name_bytes, count, atom_per_mol
                    )
                )

            self.__write(struct.pack("<I", self.n_atoms))
            for name in res_names:
                name_bytes = name.encode("utf-8")
                self.__write(
                    struct.pack(f"<B{len(name_bytes)}s", len(name_bytes), name_bytes)
                )

            for res_id in res_ids:
                self.__write(struct.pack("<I", res_id))
            self.__flush()
        else:
            # verify consistent magic number
            self.__handle = open(path, "rb+")  # noqa: SIM115
            magic_number = self.__handle.read(8)
            if magic_number != magic_1_0:
                raise ValueError(
                    "Magic number mismatch. Can't append to a different file format."
                )

            # truncate if needed
            if truncate is not None:
                # skip header
                (header_size,) = struct.unpack("<Q", self.__handle.read(8))
                self.__handle.seek(8 + header_size + 4, SEEK_CUR)

                # start reading frames
                truncate_at = self.__handle.tell()
                while chunk_header := self.__handle.read(8):
                    (chunk_size,) = struct.unpack("<Q", chunk_header)
                    self.__handle.seek(8, SEEK_CUR)
                    content = zlib.decompress(self.__handle.read(chunk_size))
                    frame, n_atoms, sim_step, sim_time = struct.unpack(
                        "<IIQd", content[0:24]
                    )
                    assert frame == self.__last_frame + 1
                    self.__handle.seek(4, SEEK_CUR)
                    if sim_step > truncate:
                        break
                    truncate_at = self.__handle.tell()
                    self.__last_frame += 1
                self.__handle.truncate(truncate_at)

            # reopen for appending
            self.__handle.close()
            self.__handle = open(path, "ab")

    def __enter__(self) -> TopTrajWriter:
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback) -> None:
        self.finish()

    def new_frame(
        self, frame_num: int, sim_step: int, time_ps: float, n_atoms: int
    ) -> None:
        """Write the frame header to the disk for a new frame.

        Frames should be written in this order:

        - new_frame()
        - write_frame_atoms()
        - write_frame_bonds()
        - write_frame()

        :param frame_num: The current frame number. Must be one larger than the previous frame.
        :param sim_step: The simulation step.
        :param time_ps: The simulation time in picoseconds.
        :param n_atoms: The number of atoms in this frame.
        """
        assert frame_num - self.__last_frame == 1, (
            "Frames passed to TopTrajWriter must be in a sequence. "
            f"Got frame num {frame_num}, expected {self.__last_frame + 1}."
        )
        assert self.__frame is None, "Only call new_frame after write_frame()!"
        self.__last_frame = frame_num

        self.__write(struct.pack("<IIQd", frame_num, n_atoms, sim_step, time_ps))

        self.__frame = {"n_atoms": n_atoms}

    def write_frame_atoms(
        self,
        names: Collection[str],
        atom_types: Collection[str],
        charges: Collection[float],
        masses: Collection[float],
    ) -> None:
        """Writes the current frame atom information to disk.

        The number of atoms must be the same as n_atoms specified in new_frame().
        new_frame() must be called first. register_frame_atoms() must be called exactly once per frame.

        :param names: The names of the atoms.
        :param atom_types: The Non-Bonded force atom types, per atom.
        :param charges: The charges of the atoms, per atom.
        :param masses: The masses of the atoms, per atom.
        """
        if self.__frame is None:
            raise ValueError("Must call new_frame() first!")

        if self.__frame.get("atoms") is not None:
            raise ValueError("Must only call write_frame_atoms once per frame!")

        if self.__frame.get("bonds") is not None:
            raise ValueError(
                "Must only call write_frame_atoms before write_frame_bonds()!"
            )

        assert (
            len(names) == len(atom_types) == len(charges) == len(masses) == self.n_atoms
        )
        assert len(names) == self.__frame["n_atoms"]

        for name in names:
            name_bytes = name.encode("utf-8")
            self.__write(
                struct.pack(f"<B{len(name_bytes)}s", len(name_bytes), name_bytes)
            )

        for atom_type in atom_types:
            atom_type_bytes = atom_type.encode("utf-8")
            self.__write(
                struct.pack(
                    f"<B{len(atom_type_bytes)}s", len(atom_type_bytes), atom_type_bytes
                )
            )

        for charge in charges:
            self.__write(struct.pack("<f", charge))

        for mass in masses:
            self.__write(struct.pack("<f", mass))

        self.__frame["atoms"] = True

    def write_frame_bonds(self, bonds: BondGraph) -> None:
        """Writes the bonds for the current frame to disk.

        :param bonds: The bonds to write, as a BondGraph object.
        """
        if self.__frame is None:
            raise ValueError("Must call new_frame() first!")
        if self.__frame.get("bonds") is not None:
            raise ValueError("Must only call write_frame_bonds once per frame!")
        if self.__frame.get("atoms") is None:
            raise ValueError("Must call write_frame_atoms before write_frame_bonds()!")

        self.__frame["bonds"] = True
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
        self.__write(raw_bytes)

    def write_frame(self) -> None:
        """Finishes writing the current frame to disk.

        Must call register_frame_atoms and register_frame_bonds exactly once first.

        Note: will automatically flush after finishing the frame.
        """
        self.__frame = None
        self.__flush()

    def __write(self, raw_bytes: bytes | bytearray) -> None:
        """Writes raw bytes to internal buffer."""
        self.__buffer += raw_bytes

    def __flush(self) -> None:
        """Writes the internal buffer to disk.

        * Uses zlib for compression.
        * Prepends the number of compressed bytes.
        * Appends the CRC32 of the compressed contents.
        """
        assert self.__handle is not None
        assert len(self.__buffer) > 0
        crc32 = zlib.crc32(self.__buffer)
        compressed = zlib.compress(self.__buffer)
        self.__handle.write(struct.pack("<Q", len(compressed)))
        self.__handle.write(struct.pack("<Q", len(self.__buffer)))
        self.__handle.write(compressed)
        self.__handle.write(struct.pack("<I", crc32))
        self.__handle.flush()
        self.__buffer = bytearray()

    def finish(self) -> None:
        """Must call when done writing the file."""
        assert self.__handle is not None
        assert len(self.__buffer) == 0, "finish() called when a frame is not finished."
        self.__handle.close()
        self.__handle = None


@dataclass
class TopTrajFrame:
    frame_index: int
    sim_step: int
    sim_time: float
    n_atoms: int
    names: list[bytes]
    res_names: list[bytes]
    res_ids: list[int]
    atom_types: list[bytes]
    charges: list[float]
    masses: list[float]
    bonds: list[tuple[int, int]]


class TopTrajReader:
    def __init__(self, path: str) -> None:
        """Create a ``.toptraj`` file reader."""
        self.path = path
        self.__handle = open(path, "rb")
        self.header = self.__handle.read(8)
        assert self.header == b"\xc0TOPTR\x01\x00"
        self.major_version = 1
        self.minor_version = 0

        header = self.__read_chunk()
        if header is None:
            raise ValueError(
                "No header was found in .toptraj file, therefore it is likely corrupt."
            )

        title_len = header[0]
        (self.title,) = struct.unpack(f"<{title_len}s", header[1 : 1 + title_len])
        i = 1 + title_len
        (n_init,) = struct.unpack("<I", header[i : i + 4])
        i += 4
        self.initial_molecules: list[tuple[bytes, int, int]] = []
        for _ in range(n_init):
            name_len = header[i]
            i += 1
            (name,) = struct.unpack(f"<{name_len}s", header[i : i + name_len])
            i += name_len
            (n, n_atoms_per) = struct.unpack("<II", header[i : i + 8])
            i += 8
            self.initial_molecules.append((name, n, n_atoms_per))

        (n_atoms,) = struct.unpack("<I", header[i : i + 4])
        i += 4
        i, self.res_names = self.__read_strings(header, i, n_atoms)
        i, self.res_ids = self.__read_any(header, i, n_atoms, "<I", 4)
        assert i == len(header)
        del header
        self.frame = 0

    def __enter__(self) -> TopTrajReader:
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback) -> None:
        self.close()

    def skip_frame(self) -> None:
        """Skip a single frame.

        Note: if already at the end, will silently do nothing.
        """
        chunk_len_bytes = self.__handle.read(8)
        if len(chunk_len_bytes) == 0:
            return
        (chunk_len,) = struct.unpack("<Q", chunk_len_bytes)
        self.__handle.seek(8 + chunk_len + 4, SEEK_CUR)

    def tell(self) -> int:
        """Return the current position in the file."""
        return self.__handle.tell()

    def __read_chunk(self) -> bytes | None:
        chunk_len_bytes = self.__handle.read(8)
        if len(chunk_len_bytes) == 0:
            # EOF
            return None
        assert len(chunk_len_bytes) == 8, (
            f"File {self.path} is likely corrupt. Trailing bytes found."
        )
        (chunk_len,) = struct.unpack("<Q", chunk_len_bytes)
        (decompressed_len,) = struct.unpack("<Q", self.__handle.read(8))
        content = zlib.decompress(self.__handle.read(chunk_len))
        assert len(content) == decompressed_len, (
            f"File {self.path} appears to be corrupt. "
            + f"Decompressed length {len(content)} doesn't match {decompressed_len}."
        )
        (crc32,) = struct.unpack("<I", self.__handle.read(4))
        if zlib.crc32(content) != crc32:
            raise ValueError(
                f"File {self.path} appears to be corrupt. Chunk CRC32 {crc32} doesn't match {zlib.crc32(content)}."
            )
        return content

    def __read_strings(
        self, content: bytes, i: int, n_atoms: int
    ) -> tuple[int, list[bytes]]:
        """Reads n pascal strings."""
        res = []
        for _ in range(n_atoms):
            len_ = content[i]
            i += 1
            res.append(struct.unpack(f"<{len_}s", content[i : i + len_])[0])
            i += len_
        return i, res

    def __read_any(
        self, content: bytes, i: int, n_atoms: int, byte_format: str, n_bytes: int
    ) -> tuple[int, list[Any]]:
        """Reads n occurrences of the byte format specified."""
        res = []
        for _ in range(n_atoms):
            res.append(struct.unpack(byte_format, content[i : i + n_bytes])[0])
            i += n_bytes
        return i, res

    def read_frame(self) -> TopTrajFrame | None:
        """Reads a frame from the ``toptraj`` file.

        :return: A toptraj frame, or None if finished.
        """
        content = self.__read_chunk()
        if content is None:
            return None

        i = 0
        frame, n_atoms, sim_step, sim_time = struct.unpack("<IIQd", content[i : i + 24])
        i += 24
        i, names = self.__read_strings(content, i, n_atoms)
        i, types = self.__read_strings(content, i, n_atoms)
        i, charges = self.__read_any(content, i, n_atoms, "<f", 4)
        i, masses = self.__read_any(content, i, n_atoms, "<f", 4)
        (n_bonds,) = struct.unpack("<Q", content[i : i + 8])
        bonds = []
        i += 8
        for _ in range(n_bonds):
            bonds.append(struct.unpack("<II", content[i : i + 8]))
            i += 8

        assert i == len(content)
        del content

        return TopTrajFrame(
            frame,
            sim_step,
            sim_time,
            n_atoms,
            names,
            self.res_names,
            self.res_ids,
            types,
            charges,
            masses,
            bonds,
        )

    def close(self) -> None:
        self.__handle.close()
