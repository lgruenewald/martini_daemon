"""Test the in-house .toptraj format reader and writer bundled with Martini Daemon."""

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

    print("generating data")
    # generate data
    n_frames = 100
    frames = [{} for _ in range(n_frames)]

    n_atoms = 5000
    res_names = [get_random_name() for _ in range(n_atoms)]
    res_ids = list(np.random.randint(1, high=n_atoms, size=n_atoms))

    for frame_num, frame in enumerate(frames):
        frame["n_atoms"] = n_atoms
        frame["names"] = [get_random_name() for _ in range(n_atoms)]
        frame["types"] = [get_random_name() for _ in range(n_atoms)]
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
    tmp_file = ".out.toptraj"
    with TopTrajWriter(
        tmp_file,
        "example title",
        [("a", 1, 1), ("b", 2, 1), ("c", n_atoms - 3, 1)],
        res_names,
        res_ids,
    ) as w:
        for frame_num, frame in enumerate(frames):
            w.new_frame(
                frame_num,
                frame_num * 5000,
                frame_num * 0.1,
                frame["n_atoms"],
            )
            w.write_frame_atoms(
                frame["names"],
                frame["types"],
                frame["charges"],
                frame["masses"],
            )
            w.write_frame_bonds(frame["bonds"])
            w.write_frame()

    print("reading from file and verifying")
    with TopTrajReader(tmp_file) as r:
        assert r.title == "example title"
        assert r.initial_molecules[0] == ("a", 1, 1)
        assert r.initial_molecules[1] == ("b", 2, 1)
        assert r.initial_molecules[2] == ("c", n_atoms - 3, 1)

        for j in range(r.n_atoms):
            assert r.res_names[j] == res_names[j]
            assert r.res_ids[j] == res_ids[j]

        for i, frame in enumerate(frames):
            if i % 5 == 0:
                # also test skipping frames
                r.skip_frame()
                continue
            f = r.read_frame()
            assert f is not None
            assert f.sim_step == i * 5000
            assert f.frame_index == i
            assert isclose(f.sim_time, i * 0.1, rel_tol=1e-5)
            assert f.n_atoms == frame["n_atoms"]
            for j in range(f.n_atoms):
                assert f.names[j] == frame["names"][j]
                assert f.atom_types[j] == frame["types"][j]
                assert isclose(f.charges[j], frame["charges"][j], rel_tol=1e-5)
                assert isclose(f.masses[j], frame["masses"][j], rel_tol=1e-5)
            assert len(f.bonds) == len(frame["bonds"].to_list())
            # Note: if this fails in the future, consider if the order (i, j) and (j, i) is not inverted
            # currently there is a fixed order in to_list, which the writer also calls, so a simple equality is fine
            assert len(set(f.bonds) - set(frame["bonds"].to_list())) == 0
            assert len(set(frame["bonds"].to_list()) - set(f.bonds)) == 0

        assert r.read_frame() is None

    os.remove(tmp_file)
