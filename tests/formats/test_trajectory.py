import os

import numpy as np
import pytest

from martini_daemon import PeriodicBox, TrajectoryReader, TrajectoryWriter


@pytest.fixture
def rootdir(request):
    return os.path.dirname(request.path)


backends = TrajectoryWriter.backends


@pytest.mark.parametrize("backend", backends)
def test_read_write_trajectory(backend: str, rootdir: str) -> None:
    """Read/write trajectories in available backends.

    Will skip formats with missing optional dependencies.
    """
    os.chdir(rootdir)
    format = "." + backend.split("_")[0]
    tmp_path = ".tmp" + format
    n_frames = 1000
    n_atoms = 2000

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
        w.write_frame(i * 100000000, i * 5.0, pbc, pos[i], vel[i])

    w.finish()

    r = TrajectoryReader(tmp_path, backend)
    for i in range(n_frames):
        f = r.read_frame()
        assert f is not None
        (sim_step, sim_time, sim_pbc, sim_pos, sim_vel) = f
        assert sim_step == i * 100000000
        assert np.isclose(sim_time, 5.0 * i, atol=1e-3)
        assert np.allclose(sim_pbc.a, pbc.a, atol=1e-3)
        assert np.allclose(sim_pbc.b, pbc.b, atol=1e-3)
        assert np.allclose(sim_pbc.c, pbc.c, atol=1e-3)
        assert np.allclose(pos[i], sim_pos, atol=1e-3)
        if sim_vel is not None:
            assert np.allclose(vel[i], sim_vel, atol=1e-3)

    os.remove(tmp_path)
