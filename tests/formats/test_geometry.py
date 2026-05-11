"""Test the different geometry file readers and writers bundled with Martini Daemon."""

import os

import numpy as np
import pytest

from martini_daemon import PeriodicBox, read_geometry, write_geometry


@pytest.fixture
def rootdir(request: pytest.FixtureRequest) -> str:
    """Fixture to get root directory."""
    return os.path.dirname(request.path)


formats = [".xyz", ".gro"]


@pytest.mark.parametrize("format_", formats)
def test_read_write_geometry(format_: str, rootdir: str) -> None:
    """Test by writing random data to a file and reading it back."""
    os.chdir(rootdir)
    tmp_path = ".tmp" + format_
    # intentionally more than 5 characters
    n_atoms = 200000

    names = ["W"] * n_atoms
    resnames = ["W"] * n_atoms
    resids = list(range(n_atoms))

    box_size = 10.0
    pbc = PeriodicBox.cubic(box_size)
    pos = np.random.rand(n_atoms, 3) * box_size
    max_vel = 50.0
    vel = np.random.rand(n_atoms, 3) * max_vel * 2.0 - max_vel

    write_geometry(
        tmp_path,
        "this is a temporary file, you may delete it if it exists",
        names,
        resnames,
        resids,
        pbc,
        pos,
        vel,
    )

    read_box, read_pos, read_vel = read_geometry(tmp_path)

    assert read_vel is not None

    pos = pos.flatten()
    read_pos = read_pos.flatten()
    vel = vel.flatten()
    read_vel = read_vel.flatten()

    # .gro doesn't contain that many significant digits to be fair
    assert np.allclose(
        [pbc.a, pbc.b, pbc.c], [read_box.a, read_box.b, read_box.c], rtol=1e-5
    )
    largest = np.argmax(np.abs(read_pos - pos))
    assert np.allclose(read_pos, pos, atol=1e-3), (
        f"Largest deviation in pos at {largest // 3}, with a difference of {np.abs(read_pos[largest] - pos[largest])} between original {pos[largest]} and read back {read_pos[largest]}."
    )
    largest = np.argmax(np.abs(read_vel - vel))
    assert np.allclose(read_vel, vel, atol=1e-3), (
        f"Largest deviation in pos at {largest // 3}, , with a difference of {np.abs(read_vel[largest] - vel[largest])} between original {vel[largest]} and read back {read_vel[largest]}."
    )

    os.remove(tmp_path)
