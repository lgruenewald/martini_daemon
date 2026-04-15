"""Test the in-house .toptraj format reader and writer bundled with Martini Daemon."""

#!/usr/bin/env python3

import os
import random
import string
from math import isclose

import numpy as np
import pytest

from martini_daemon import BondGraph, TopTrajReader, TopTrajWriter


def get_random_name() -> str:
    """Get a random name for atom names, residue names or atom types."""
    return random.choice(list(string.ascii_uppercase))


@pytest.fixture
def rootdir(request: pytest.FixtureRequest) -> str:
    """Fixture to get root directory."""
    return os.path.dirname(request.path)


def test_toptraj_writer(rootdir: str) -> None:
    """Write a .toptraj file with random data and read it back to compare."""
    os.chdir(rootdir)
    tmp_file = ".out.toptraj"
    w = TopTrajWriter(tmp_file, "example title", [("a", 1), ("b", 2), ("c", 3)])

    print("generating data")
    # generate data
    n_frames = 100
    frames = [{} for _ in range(n_frames)]

    for frame_num, frame in enumerate(frames):
        n_atoms = np.random.randint(2000, 3000)
        frame["n_atoms"] = n_atoms
        frame["names"] = [get_random_name() for _ in range(n_atoms)]
        frame["types"] = [get_random_name() for _ in range(n_atoms)]
        frame["res_names"] = [get_random_name() for _ in range(n_atoms)]
        frame["res_ids"] = np.random.randint(1, high=n_atoms, size=n_atoms)
        frame["charges"] = np.random.rand(n_atoms) * 2.0 - 1.0
        frame["masses"] = np.random.rand(n_atoms) * 50.0
        n_bonds = np.random.randint(1000, 5000)
        bonds = BondGraph(n_atoms)
        for _ in range(n_bonds):
            i = np.random.randint(0, n_atoms)
            j = np.random.randint(0, n_atoms)
            # BondGraph ignores it if i==j
            bonds.add_bond(i, j)
        frame["bonds"] = bonds
        # and write it to file

    print("writing it to file")
    for frame_num, frame in enumerate(frames):
        w.new_frame(
            frame_num,
            frame_num * 5000,
            frame_num * 0.1,
            frame["n_atoms"],
        )
        w.write_frame_atoms(
            frame["names"],
            frame["res_names"],
            frame["res_ids"],
            frame["types"],
            frame["charges"],
            frame["masses"],
        )
        w.write_frame_bonds(frame["bonds"])
        w.write_frame()
    w.finish()

    print("reading from file and verifying")
    r = TopTrajReader(tmp_file)
    assert r.title == b"example title"
    assert r.initial_molecules[0] == (b"a", 1)
    assert r.initial_molecules[1] == (b"b", 2)
    assert r.initial_molecules[2] == (b"c", 3)

    for i, frame in enumerate(frames):
        f = r.read_frame()
        assert f is not None
        assert f.sim_step == i * 5000
        assert f.frame_index == i
        assert isclose(f.sim_time, i * 0.1, rel_tol=1e-5)
        assert f.n_atoms == frame["n_atoms"]
        for j in range(f.n_atoms):
            assert f.names[j] == frame["names"][j].encode("utf-8")
            assert f.res_names[j] == frame["res_names"][j].encode("utf-8")
            assert f.atom_types[j] == frame["types"][j].encode("utf-8")
            assert isclose(f.charges[j], frame["charges"][j], rel_tol=1e-5)
            assert isclose(f.masses[j], frame["masses"][j], rel_tol=1e-5)
            assert f.res_ids[j] == frame["res_ids"][j]
        assert len(f.bonds) == len(frame["bonds"].to_list())
        # Note: if this fails in the future, consider if the order (i, j) and (j, i) is not inverted
        # currently there is a fixed order in to_list, which the writer also calls, so a simple equality is fine
        assert len(set(f.bonds) - set(frame["bonds"].to_list())) == 0
        assert len(set(frame["bonds"].to_list()) - set(f.bonds)) == 0

    assert r.read_frame() is None

    os.remove(tmp_file)
