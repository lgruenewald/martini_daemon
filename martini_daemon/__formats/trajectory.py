from os.path import splitext

import numpy as np
import numpy.typing as npt
from openmm.app.internal.xtc_utils import (  # ty: ignore[unresolved-import]
    read_xtc,
    xtc_write_frame,
)

from ..__rust import PeriodicBox


def _to_int32(n: int) -> int:
    """Put a value back into int32"""
    n = n % (np.iinfo(np.uint32).max + 1)
    if n > np.iinfo(np.int32).max:
        n -= np.iinfo(np.uint32).max + 1
    return n


def _from_int32(n: int, last: int) -> int:
    n = int(n)
    if n < 0:
        n += np.iinfo(np.uint32).max + 1
    while n < last:
        n += np.iinfo(np.uint32).max + 1
    return n


class TrajectoryWriter:
    backends = [
        "xtc_openmm_internal",
        "xtc_molly",
        "trr_mdtraj",
    ]

    default_backends = {
        ".xtc": "xtc_openmm_internal",
        ".trr": "trr_mdtraj",
    }

    def __init__(
        self,
        path: str,
        backend: str | None = None,
        append: bool = False,
        truncate: int | None = None,
    ) -> None:
        """Create a TrajectoryWriter object.

        Note: based on the backend, this object may own open file handles.
        These must be closed manually with .finish() if using this class!

        :param path: Path to the trajectory file.
        :param backend: Trajectory writer backend. If None, writer will be determined based on the file extension.
            See TrajectoryWriter.backends for the list of available backends, as well as TrajectoryWriter.default_backends
            for the default choices. Some backends may require self-explanatory optional dependencies.
            See pyproject.toml or the README in the repo for details.
        :param append: Whether to append to the trajectory or not.
        :param truncate: If appending, the last simulation step to keep.
        """
        self.path = path
        if not append:
            assert truncate is None, "Truncate is only valid if append is True."
        if backend is None:
            _, ext = splitext(path)
            self.backend = self.default_backends.get(ext)
            if self.backend is None:
                raise ValueError(f"Unknown trajectory file format {ext} for {path}.")
        else:
            self.backend = backend
        if self.backend not in self.backends:
            raise ValueError(
                f"Unknown trajectory backend {self.backend}. Available Trajectory Writer backends: {', '.join(self.backends)}."
            )
        match self.backend:
            case "xtc_openmm_internal":
                # reopens the file every time I guess
                assert append is False, "can't append with xtc_openmm_internal"
                assert truncate is None, (
                    "Can't truncate file using xtc_openmm_internal."
                )
                pass
            case "xtc_molly":
                import molly

                if append and truncate is not None:
                    reader = molly.XTCReader(path)
                    last = 0
                    last_tell = 0
                    while last < truncate:
                        last_tell = reader.tell()
                        f = reader.pop_frame()
                        last = _from_int32(f.step, last)
                    reader.close()

                    with open(path, "rb+") as f:
                        # f.seek(last_tell)
                        f.truncate(last_tell)
                        # assert f.tell() == last_tell

                self.__writer_molly = molly.XTCWriter(path, append)
            case "trr_mdtraj":
                import mdtraj.formats

                assert truncate is None, "TODO"  # TODO
                assert append is False, "TODO"  # TODO

                self.__writer_trr = mdtraj.formats.TRRTrajectoryFile(path, "w")
            case _:
                assert False

    def write_frame(
        self,
        sim_step: int,
        time_ps: float,
        box: PeriodicBox,
        pos: npt.NDArray[np.float32 | np.float64],
        vel: npt.NDArray[np.float32 | np.float64] | None = None,
    ) -> None:
        """Write a single frame to the trajectory file.

        :param sim_step: Current simulation step.
        :param time_ps: Current simulation time in ps.
        :param box: Current box in nanometers.
        :param pos: Current position in nanometers.
        :param vel: Current velocity in nanometers / picosecond. Will be only written if the format supports it.
        """
        if vel is not None:
            assert vel.shape == pos.shape

        match self.backend:
            case "xtc_openmm_internal":
                pos = np.array(pos, dtype=np.float32)
                n_atoms = len(pos)
                assert pos.shape == (n_atoms, 3)
                box_numpy = np.array([box.a, box.b, box.c], dtype=np.float32)
                assert box_numpy.shape == (3, 3)
                if pos.dtype != np.float32:
                    pos = np.array(pos, dtype=np.float32)
                xtc_write_frame(
                    self.path.encode("utf-8"),  # title as byte string
                    pos,  # positions as float[:, :]
                    box_numpy,  # box as float[:, :]
                    time_ps,  # time in ps
                    # _to_int32(sim_step),
                    sim_step,
                )
            case "xtc_molly":
                import molly

                frame = molly.Frame()
                frame.box = [box.a, box.b, box.c]
                frame.positions = pos.flatten()
                frame.step = sim_step % (np.iinfo(np.uint32).max + 1)
                frame.time = time_ps
                assert self.__writer_molly is not None
                self.__writer_molly.write_frame(frame)
            case "trr_mdtraj":
                assert len(pos.shape) == 2 and pos.shape[1] == 3
                self.__writer_trr._write(
                    np.array([pos], dtype=np.float32),
                    np.array([time_ps], dtype=np.float32),
                    np.array([_to_int32(sim_step)], dtype=np.int32),
                    np.array([[box.a, box.b, box.c]], dtype=np.float32),
                    np.zeros(1, dtype=np.float32),
                    np.array([vel], dtype=np.float32) if vel is not None else None,
                    None,
                )
            case _:
                assert False

    def finish(self) -> None:
        """Close any open file handles, depending on the backend."""
        match self.backend:
            case "xtc_openmm_internal":
                pass
            case "xtc_molly":
                assert self.__writer_molly is not None
                self.__writer_molly.close()
                self.__writer_molly = None
            case "trr_mdtraj":
                self.__writer_trr.close()
            case _:
                assert False


class TrajectoryReader:
    backends = [
        "xtc_openmm_internal",
        "xtc_molly",
        "trr_mdtraj",
    ]

    default_backends = {
        ".xtc": "xtc_openmm_internal",
        ".trr": "trr_mdtraj",
    }

    default_backends = {".xtc": "xtc_openmm_internal"}

    def __init__(self, path: str, backend: str | None = None) -> None:
        """Create a TrajectoryReader object.

        Note: based on the backend, this object may own open file handles.
        These must be closed manually with .finish() if using this class!

        :param path: Path to the trajectory file.
        :param backend: Trajectory reader backend. If None, writer will be determined based on the file extension.
            Current list of backends: "xtc_openmm_internal".
        """
        self.path = path
        if backend is None:
            _, ext = splitext(path)
            self.backend = self.default_backends.get(ext)
            if self.backend is None:
                raise ValueError(f"Unknown trajectory file format {ext} for {path}.")
        if backend not in self.backends:
            raise ValueError(
                f"Unknown trajectory backend {backend}. Available Trajectory Writer backends: {', '.join(self.backends)}."
            )
        self.backend = backend
        match self.backend:
            case "xtc_openmm_internal":
                self.__pos, self.__box, self.__time, self.__step = read_xtc(
                    path.encode("utf-8")  # must be bytestring
                )
                self.__n_frames = self.__pos.shape[2]
                self.__c_frame = 0
                self.__last_step = 0
            case "xtc_molly":
                import molly

                self.__reader_molly = molly.XTCReader(path)
                self.__last_step = 0
            case "trr_mdtraj":
                import mdtraj.formats

                self.__last_step = 0
                self.__reader_trr = mdtraj.formats.TRRTrajectoryFile(path, "r")
            case _:
                assert False

    def read_frame(
        self,
    ) -> None | tuple[int, float, PeriodicBox, np.ndarray, np.ndarray | None]:
        """Read a single frame from the trajectory file.

        :return: Current simulation step, simulation time in ps, periodic box, positions in nm
            and the velocities in nm/ps if format supports it, or None if not. None if all frames were read.
        """
        match self.backend:
            case "xtc_openmm_internal":
                if self.__c_frame >= self.__n_frames:
                    return None
                self.__c_frame += 1
                step = _from_int32(self.__step[self.__c_frame - 1], self.__last_step)
                self.__last_step = step

                return (
                    step,
                    self.__time[self.__c_frame - 1],
                    PeriodicBox(
                        self.__box[:, 0, self.__c_frame - 1],
                        self.__box[:, 1, self.__c_frame - 1],
                        self.__box[:, 2, self.__c_frame - 1],
                    ),
                    self.__pos[:, :, self.__c_frame - 1],
                    None,
                )
            case "xtc_molly":
                assert self.__reader_molly is not None
                self.__reader_molly.read_frame()
                if self.__reader_molly.frame is None:
                    return None
                step = _from_int32(self.__reader_molly.frame.step, self.__last_step)
                self.__last_step = step
                return (
                    step,
                    self.__reader_molly.frame.time,
                    PeriodicBox(*self.__reader_molly.frame.box),
                    self.__reader_molly.frame.positions,
                    None,
                )
            case "trr_mdtraj":
                assert self.__reader_trr is not None
                try:
                    xyz, time, step, box, lambd, vel, force = self.__reader_trr._read(
                        1, None, True
                    )
                    assert vel is not None
                except RuntimeError:
                    # no velocities in file
                    xyz, time, step, box, lambd, vel, force = self.__reader_trr._read(
                        1, None
                    )
                    assert vel is None
                if len(xyz) == 0:
                    return None
                step = _from_int32(step[0], self.__last_step)
                self.__last_step = step
                return (
                    step,
                    time,
                    PeriodicBox(box[0, 0], box[0, 1], box[0, 2]),
                    xyz[0],
                    vel[0] if vel is not None else None,
                )

            case _:
                assert False

    def finish(self) -> None:
        """Close any open file handles, depending on the backend."""
        match self.backend:
            case "xtc_openmm_internal":
                pass
            case "xtc_molly":
                self.__reader_molly.close()
            case "trr_mdtraj":
                self.__reader_trr.close()
            case _:
                assert False
