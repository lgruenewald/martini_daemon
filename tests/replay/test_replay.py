"""Replay algorithm tests."""

import glob
import os

import numpy as np
import pytest

from martini_daemon import (
    CheckpointLoader,
    CheckpointReporter,
    PeriodicBox,
    ReactionReporter,
    Simulation,
    SystemDump,
)

e_tol = 1e-5  # energy relative tolerance
f_tol = 2e-4  # force relative tolerance


# == CONFIG ==
@pytest.fixture
def rootdir(request: pytest.FixtureRequest) -> str:
    """Get the root directory for this test."""
    return os.path.dirname(request.path)


tests = [
    "update_redefine",
]


def load_checkpoint(
    checkpoint: str, top: str, gro: str
) -> tuple[float, np.ndarray, PeriodicBox, np.ndarray, np.ndarray]:
    """Load a checkpoint and return the energy and forces"""
    sim = CheckpointLoader(
        checkpoint,
        top,
        gro,
        0,
        platform="Reference",
    )
    # manually dump
    d = SystemDump()
    d.on_simulation_start(sim, False)
    d.on_simulation_finish(sim)

    assert sim.current_step == 1
    assert sim.reactions_so_far == 1
    ke, pe, te = sim.context.get_energies()
    forces = sim.context.get_forces()
    pos, box = sim.context.get_positions()
    vel = sim.context.get_velocities()
    sim.finish()
    return pe, forces, box, pos, vel


def get_chk(top: str, gro: str) -> tuple[str, PeriodicBox, np.ndarray, np.ndarray]:
    """Get a simulation checkpoint based on a single modification frame, based on top and gro."""
    sim = Simulation(top, gro, 0, reporters=[CheckpointReporter(), ReactionReporter()])
    sim.step(0, traj=False, dm=True)
    sim.current_step = 1
    sim.step(0, traj=True, dm=False)
    pos, box = sim.context.get_positions()
    vel = sim.context.get_velocities()
    sim.finish()
    return "out.chk", box, pos, vel


def run_gromacs(x: str) -> tuple[float, np.ndarray]:
    """Get GROMACS energy and force for the system in the current directory."""
    os.system("../gmxrun.sh")
    assert os.path.isfile("energy.xvg"), f"./gmxrun.sh failure for {x} (E)"
    assert os.path.isfile("forces.xvg"), f"./gmxrun.sh failure for {x} (F)"
    with open("energy.xvg") as f:
        lines = f.readlines()
        gmx_energy = float(lines[-1].split()[-1])
    with open("forces.xvg") as f:
        lines = f.readlines()
        gmx_force_line = lines[-1].split()
        gmx_forces = np.array([float(x) for x in gmx_force_line][1:])
    return gmx_energy, gmx_forces


@pytest.mark.parametrize("x", tests)
def test_replay(x: str, rootdir: str) -> None:
    """Test the modification algorithm by comparing SysStar dumps."""
    # enter dir
    os.chdir(rootdir)
    assert os.path.isdir(x)
    assert os.path.isfile("gmxrun.sh")
    os.chdir(x)

    chk, old_box, old_pos, old_vel = get_chk("system.top", "system.gro")
    energy, forces, box, pos, vel = load_checkpoint(chk, "system.top", "system.gro")
    assert np.allclose(old_box.a, box.a)
    assert np.allclose(old_box.b, box.b)
    assert np.allclose(old_box.c, box.c)
    assert np.allclose(old_pos, pos)
    assert np.allclose(old_vel, vel)
    gmx_energy, gmx_forces = run_gromacs(x)
    assert np.isclose(energy, gmx_energy, rtol=e_tol, atol=0)
    a_diffs = np.abs(forces.flatten() - gmx_forces)
    r_diffs = a_diffs / np.abs(gmx_forces)
    largest_diff = np.argmax(r_diffs)
    assert np.all(r_diffs < f_tol), (
        f"Force deviation maximum on atom {largest_diff // 3} (0 indexed), dimension {largest_diff % 3}\n"
        + f" abs diff {a_diffs[largest_diff]}\n"
        + f" rel diff {r_diffs[largest_diff]}\n"
        + f" tolerance: {(f_tol * np.abs(gmx_forces[largest_diff]))}\n"
        + f" gmx: {gmx_forces[largest_diff]}\n martini daemon: {forces.flatten()[largest_diff]}."
    )

    # cleanup
    for filename in glob.glob("./out*"):
        os.remove(filename)
    for filename in glob.glob("./#*"):
        os.remove(filename)
    os.remove("energy.xvg")
    os.remove("forces.xvg")
