"""
Test the daemon cli's operations.

Note: This testing only runs them to see if the exit code is 0.
"""

import os
from subprocess import run

import pytest

from martini_daemon import TopTrajReader, TrajectoryReader


@pytest.fixture
def rootdir(request: pytest.FixtureRequest) -> str:
    """Get the root directory where this test is located."""
    return os.path.dirname(request.path)


def test_info(rootdir: str) -> None:
    """Test the daemon info command."""
    os.chdir(rootdir)
    run(["daemon", "info", "out.chk"], check=True)
    run(["daemon", "info", "sel2.toptraj"], check=True)


def test_select(rootdir: str) -> None:
    """Test the daemon select command on toptraj."""
    os.chdir(rootdir)
    run(
        [
            "daemon",
            "select",
            "-i",
            "sel2.toptraj",
            "-o",
            "tmp.toptraj",
            "-f",
            "::5",
            "-a",
            ":100",
        ],
        check=True,
    )
    os.unlink("tmp.toptraj")


def test_whole(rootdir: str) -> None:
    """Test the daemon whole command using a toptraj file."""
    os.chdir(rootdir)
    run(
        ["daemon", "whole", "-i", "sel2.xtc", "-o", "tmp.xtc", "-s", "sel2.toptraj"],
        check=True,
    )
    with TrajectoryReader("tmp.xtc") as xtc, TopTrajReader("sel2.toptraj") as top:
        while (f := xtc.read_frame()) is not None:
            tf = top.read_frame()
            assert tf is not None
            _, _, box, pos, _ = f
            for i, j in tf.bonds:
                assert not box.crosses_box(pos[i], pos[j])
    os.unlink("tmp.xtc")
