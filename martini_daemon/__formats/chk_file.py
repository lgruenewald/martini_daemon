import struct
from dataclasses import dataclass
from zlib import crc32

import numpy as np
import numpy.typing as npt

from ..__rust import PeriodicBox

# the current magic number
# increment the last byte each time the checkpoint format changes
# also increment if any reporter's output changes that we depend on -- notably ReactionReporter at the moment
# also increment if .top / .rx reaction template syntax changes
MAGIC = b"\xc0DMNCK\x00\x01"


def write_checkpoint(
    path: str,
    sim_name: str,
    current_step: int,
    trajectory_frame: int,
    reactions_so_far: int,
    time_ps: float,
    n_atoms: int,
    box: PeriodicBox,
    pos: npt.NDArray[np.float64],
    vel: npt.NDArray[np.float64],
) -> None:
    """Write a checkpoint file based on the information provided in the arguments."""
    crc = crc32(b"")
    with open(path, "wb") as f:

        def write_and_crc(b: bytes) -> None:
            nonlocal crc, f
            crc = crc32(b, crc)
            f.write(b)

        write_and_crc(MAGIC)
        sim_name_bytes = sim_name.encode("utf-8")
        write_and_crc(struct.pack("<I", len(sim_name_bytes)))
        write_and_crc(sim_name_bytes)
        write_and_crc(struct.pack("<Q", current_step))
        write_and_crc(struct.pack("<Q", trajectory_frame))
        write_and_crc(struct.pack("<Q", reactions_so_far))
        write_and_crc(struct.pack("<d", time_ps))
        write_and_crc(struct.pack("<Q", n_atoms))
        box_np = np.array([box.a, box.b, box.c], dtype=np.float64)
        assert box_np.shape == (3, 3)
        assert pos.dtype == np.float64
        assert vel.dtype == np.float64
        assert pos.shape == (n_atoms, 3)
        assert vel.shape == (n_atoms, 3)
        write_and_crc(box_np.tobytes())
        write_and_crc(pos.tobytes())
        write_and_crc(vel.tobytes())
        f.write(struct.pack("<I", crc))


@dataclass
class Checkpoint:
    sim_name: str
    current_step: int
    trajectory_frame: int
    reactions_so_far: int
    time_ps: float
    n_atoms: int
    box: PeriodicBox
    pos: npt.NDArray[np.float64]
    vel: npt.NDArray[np.float64]


def read_checkpoint(path: str) -> Checkpoint:
    with open(path, "rb") as f:
        crc = crc32(b"")

        def read_and_crc(n: int) -> bytes:
            nonlocal crc, f
            res = f.read(n)
            crc = crc32(res, crc)
            return res

        magic = read_and_crc(8)
        if magic != MAGIC:
            raise ValueError(
                f"MAGIC number mismatch, {path} is not a valid .chk file, or is from a different version."
            )
        sim_name_len = struct.unpack("<I", read_and_crc(4))[0]
        sim_name = read_and_crc(sim_name_len).decode("utf-8")

        current_step, trajectory_frame, reactions_so_far, time_ps, n_atoms = (
            struct.unpack("<QQQdQ", read_and_crc(5 * 8))
        )
        box_np = (
            np.frombuffer(read_and_crc(8 * 9), dtype=np.float64).reshape(3, 3).copy()
        )
        box = PeriodicBox(box_np[0], box_np[1], box_np[2])
        pos = (
            np.frombuffer(read_and_crc(8 * 3 * n_atoms), dtype=np.float64)
            .reshape(n_atoms, 3)
            .copy()
        )
        vel = (
            np.frombuffer(read_and_crc(8 * 3 * n_atoms), dtype=np.float64)
            .reshape(n_atoms, 3)
            .copy()
        )

        ref_crc = struct.unpack("<I", f.read(8))[0]
        if crc != ref_crc:
            raise ValueError(f"CRC mismatch, {path} is corrupt.")

        return Checkpoint(
            sim_name,
            current_step,
            trajectory_frame,
            reactions_so_far,
            time_ps,
            n_atoms,
            box,
            pos,
            vel,
        )
