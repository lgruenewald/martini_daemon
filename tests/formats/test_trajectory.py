"""Test the features of different trajectory formats supported in Martini Daemon."""

import os

import numpy as np
import pytest

from martini_daemon import PeriodicBox, TrajectoryReader, TrajectoryWriter


@pytest.fixture
def rootdir(request: pytest.FixtureRequest) -> str:
    """Return the root directory for this test."""
    return os.path.dirname(request.path)


backends = TrajectoryWriter.backends


@pytest.mark.parametrize("backend", backends)
def test_read_write_trajectory(backend: str, rootdir: str) -> None:
    """
    Read/write trajectories in available backends.

    Will skip formats with missing optional dependencies.
    """
    os.chdir(rootdir)
    format = "." + backend.split("_")[0]
    tmp_path = ".tmp" + format
    if os.path.isfile(tmp_path):
        os.remove(tmp_path)
    n_frames = 1000
    n_atoms = 2000
    # making sure it overflows int32!
    steps_per_frame = 10000000

    box_size = 10.0
    pbc = PeriodicBox.cubic(box_size)
    pos = np.random.rand(n_frames, n_atoms, 3) * box_size
    max_vel = 50.0
    vel = np.random.rand(n_frames, n_atoms, 3) * max_vel * 2.0 - max_vel

    try:
        w = TrajectoryWriter(tmp_path, backend)
    except ImportError:
        # missing optional dependency - skip test
        pytest.skip("Missing optional dependency.")
        return

    for i in range(n_frames):
        w.write_frame(i * steps_per_frame, i * 5.0, pbc, pos[i], vel[i])

    w.close()

    with TrajectoryReader(tmp_path, backend) as r:
        for i in range(n_frames):
            f = r.read_frame()
            assert f is not None
            (sim_step, sim_time, sim_pbc, sim_pos, sim_vel) = f
            assert sim_step == i * steps_per_frame
            assert np.isclose(sim_time, 5.0 * i, atol=1e-3)
            assert np.allclose(sim_pbc.a, pbc.a, atol=1e-3)
            assert np.allclose(sim_pbc.b, pbc.b, atol=1e-3)
            assert np.allclose(sim_pbc.c, pbc.c, atol=1e-3)
            assert np.allclose(pos[i], sim_pos, atol=1e-3)
            if sim_vel is not None:
                assert np.allclose(vel[i], sim_vel, atol=1e-3)
        assert r.read_frame() is None

    os.remove(tmp_path)


@pytest.mark.parametrize("backend", backends)
def test_truncate_append_trajectory(backend: str, rootdir: str) -> None:
    """
    Read/write trajectories in available backends.

    Will skip formats with missing optional dependencies.
    """
    if backend == "trr_mdtraj":
        pytest.skip("TRR mdtraj does not support appending.")
    os.chdir(rootdir)
    format = "." + backend.split("_")[0]
    tmp_path = ".tmp" + format
    if os.path.isfile(tmp_path):
        os.remove(tmp_path)

    n_frames = 1000
    truncate_at = 500
    n_atoms = 2000
    steps_per_frame = 10000000
    ps_per_frame = 5.0

    box_size = 10.0
    pbc = PeriodicBox.cubic(box_size)
    pos = np.random.rand(n_frames, n_atoms, 3) * box_size
    max_vel = 50.0
    vel = np.random.rand(n_frames, n_atoms, 3) * max_vel * 2.0 - max_vel

    try:
        w = TrajectoryWriter(tmp_path, backend)
    except ImportError:
        # missing optional dependency - skip test
        pytest.skip("Missing optional dependency.")
        return

    for i in range(n_frames):
        # making sure it overflows int32!
        w.write_frame(i * steps_per_frame, i * ps_per_frame, pbc, pos[i], vel[i])

    w.close()

    # only keep the first 500 frames
    # last frame to keep
    w = TrajectoryWriter(
        tmp_path,
        backend,
        append=True,
        keep_n_frames=truncate_at,
    )

    # rewrite the others
    pos[truncate_at:] = np.random.rand(n_frames - truncate_at, n_atoms, 3) * box_size
    vel[truncate_at:] = (
        np.random.rand(n_frames - truncate_at, n_atoms, 3) * max_vel * 2.0 - max_vel
    )

    for i in range(truncate_at, n_frames):
        w.write_frame(i * steps_per_frame, i * ps_per_frame, pbc, pos[i], vel[i])
    w.close()

    r = TrajectoryReader(tmp_path, backend)
    for i in range(n_frames):
        f = r.read_frame()
        assert f is not None
        (sim_step, sim_time, sim_pbc, sim_pos, sim_vel) = f
        assert sim_step == i * steps_per_frame
        assert np.isclose(sim_time, ps_per_frame * i, atol=1e-3)
        assert np.allclose(sim_pbc.a, pbc.a, atol=1e-3)
        assert np.allclose(sim_pbc.b, pbc.b, atol=1e-3)
        assert np.allclose(sim_pbc.c, pbc.c, atol=1e-3)
        assert np.allclose(pos[i], sim_pos, atol=1e-3)
        if sim_vel is not None:
            assert np.allclose(vel[i], sim_vel, atol=1e-3)
    assert r.read_frame() is None
    r.close()

    os.remove(tmp_path)
