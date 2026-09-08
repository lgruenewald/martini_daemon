from __future__ import annotations

import struct
import zlib
from collections.abc import Collection
from dataclasses import dataclass
from io import SEEK_CUR, SEEK_SET
from types import TracebackType
from typing import Any, Self

from ..__rust import BondGraph


class TopTrajWriter:
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
        """
        Create a Topology Trajectory writer.

        This class is the old python implementation.

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
                    frame, _n_atoms, sim_step, _sim_time = struct.unpack(
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
            self.__handle = open(path, "ab")  # noqa: SIM115

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        exc_traceback: TracebackType | None,
    ) -> None:
        self.finish()

    def new_frame(
        self, frame_num: int, sim_step: int, time_ps: float, n_atoms: int
    ) -> None:
        """
        Write the frame header to the disk for a new frame.

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
        """
        Write the current frame atom information to disk.

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

    def write_frame_bonds(self, bonds: BondGraph | list[tuple[int, int]]) -> None:
        """
        Write the bonds for the current frame to disk.

        Note: will automatically remove duplicates, self-bonds
              and will re-order bonds to have smaller first in each entry.

        :param bonds: The bonds to write, as a BondGraph object or list of (i, j) tuples.
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
        bond_list = bonds.to_list() if isinstance(bonds, BondGraph) else bonds

        for i, j in bond_list:
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
        """
        Finish writing the current frame to disk.

        Must call register_frame_atoms and register_frame_bonds exactly once first.

        Note: will automatically flush after finishing the frame.
        """
        self.__frame = None
        self.__flush()

    def __write(self, raw_bytes: bytes | bytearray) -> None:
        """Write raw bytes to internal buffer."""
        self.__buffer += raw_bytes

    def __flush(self) -> None:
        """
        Write the internal buffer to disk.

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
    names: list[str]
    atom_types: list[str]
    charges: list[float]
    masses: list[float]
    bonds: list[tuple[int, int]]


class TopTrajReader:
    def __init__(self, path: str) -> None:
        """
        Create a ``.toptraj`` file reader.

        This class is the old, python implementation.
        """
        self.path = path
        self.__handle = open(path, "rb")  # noqa: SIM115
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
        (title,) = struct.unpack(f"<{title_len}s", header[1 : 1 + title_len])
        self.title: str = title.decode("utf-8")
        i = 1 + title_len
        (n_init,) = struct.unpack("<I", header[i : i + 4])
        i += 4
        self.initial_molecules: list[tuple[str, int, int]] = []
        for _ in range(n_init):
            name_len = header[i]
            i += 1
            (name,) = struct.unpack(f"<{name_len}s", header[i : i + name_len])
            i += name_len
            (n, n_atoms_per) = struct.unpack("<II", header[i : i + 8])
            i += 8
            self.initial_molecules.append((name.decode("utf-8"), n, n_atoms_per))

        (n_atoms,) = struct.unpack("<I", header[i : i + 4])
        self.n_atoms = n_atoms
        i += 4
        i, self.res_names = self.__read_strings(header, i, n_atoms)
        i, self.res_ids = self.__read_any(header, i, n_atoms, "<I", 4)
        assert i == len(header)
        del header
        self.frame = 0

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        exc_traceback: TracebackType | None,
    ) -> None:
        self.close()

    def skip_frame(self) -> bool:
        """
        Skip a single frame.

        :return: True if a frame was skipped, False if already at the end of the file.
        """
        chunk_len_bytes = self.__handle.read(8)
        if len(chunk_len_bytes) == 0:
            return False
        (chunk_len,) = struct.unpack("<Q", chunk_len_bytes)
        self.__handle.seek(8 + chunk_len + 4, SEEK_CUR)
        return True

    def tell(self) -> int:
        """Return the current position in the file."""
        return self.__handle.tell()

    def seek(self, pos: int) -> None:
        """Set the reader to position, as returned by tell()."""
        self.__handle.seek(pos, SEEK_SET)

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
    ) -> tuple[int, list[str]]:
        """Read n pascal strings."""
        res = []
        for _ in range(n_atoms):
            len_ = content[i]
            i += 1
            bytes_: bytes
            (bytes_,) = struct.unpack(f"<{len_}s", content[i : i + len_])
            res.append(bytes_.decode("utf-8"))
            i += len_
        return i, res

    def __read_any(
        self, content: bytes, i: int, n_atoms: int, byte_format: str, n_bytes: int
    ) -> tuple[int, list[Any]]:
        """Read n occurrences of the byte format specified."""
        res = []
        for _ in range(n_atoms):
            res.append(struct.unpack(byte_format, content[i : i + n_bytes])[0])
            i += n_bytes
        return i, res

    def read_frame(self) -> TopTrajFrame | None:
        """
        Read a frame from the ``toptraj`` file.

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
            types,
            charges,
            masses,
            bonds,
        )

    def close(self) -> None:
        self.__handle.close()
