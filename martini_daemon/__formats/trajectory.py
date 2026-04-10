from os.path import splitext

import numpy as np

from ..__rust import PeriodicBox


class TrajectoryWriter:
    backends = [
        "xtc_openmm_internal",
    ]

    default_backends = {".xtc": "xtc_openmm_internal"}

    def __init__(self, path: str, backend=None):
        self.path = path
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
                pass
            case _:
                assert False

    def write_frame(
        self,
        sim_step: int,
        time_ps: float,
        box: PeriodicBox,
        pos: np.ndarray,
        vel: np.ndarray | None = None,
    ) -> None:
        match self.backend:
            case "xtc_openmm_internal":
                from openmm.app.internal.xtc_utils import (
                    xtc_write_frame,  # ty: ignore[unresolved-import]
                )

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
                    sim_step,
                )
            case _:
                assert False

    def finish(self):
        match self.backend:
            case "xtc_openmm_internal":
                pass
            case _:
                assert False


class TrajectoryReader:
    backends = [
        "xtc_openmm_internal",
    ]

    default_backends = {".xtc": "xtc_openmm_internal"}

    def __init__(self, path: str, backend=None):
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
                from openmm.app.internal.xtc_utils import (
                    read_xtc,  # ty: ignore[unresolved-import]
                )

                self.pos, self.box, self.time, self.step = read_xtc(
                    path.encode("utf-8")
                )
                self.n_frames = len(self.pos)
                self.c_frame = 0
            case _:
                assert False

    def read_frame(
        self,
    ) -> None | tuple[int, float, PeriodicBox, np.ndarray, np.ndarray | None]:
        """Returns a tuple of sim step, time (in ps), pbc, pos and vel if the format supports it"""
        match self.backend:
            case "xtc_openmm_internal":
                if self.c_frame >= self.n_frames:
                    return None
                self.c_frame += 1
                return (
                    self.step[self.c_frame - 1],
                    self.time[self.c_frame - 1],
                    PeriodicBox(
                        self.box[self.c_frame - 1, 0],
                        self.box[self.c_frame - 1, 1],
                        self.box[self.c_frame - 1, 2],
                    ),
                    self.pos[self.c_frame - 1],
                    None,
                )

            case _:
                assert False

    def finish(self):
        match self.backend:
            case "xtc_openmm_internal":
                pass
            case _:
                assert False
