#!/usr/bin/env python3

from martini_daemon import TopTrajWriter, TopTrajReader
import numpy as np
import string
from tqdm import tqdm
from math import isclose

def get_random_name():
    return "".join(
        np.random.choice(list(string.ascii_uppercase))
        for _ in range(np.random.randint(1, 8))
    )


def test_toptraj_writer():
    w = TopTrajWriter(
        "out.toptraj",
        "example title",
        [
            ("a", 1),
            ("b", 2),
            ("c", 3)
        ]
    )

    print("generating data")
    # generate data
    n_frames = 100
    frames = [
        {} for _ in range(n_frames)
    ]

    for frame_num, frame in tqdm(enumerate(frames), total=len(frames)):
        n_atoms = np.random.randint(2000, 3000)
        frame["n_atoms"] = n_atoms
        frame["names"] = [
            get_random_name() for _ in range(n_atoms)
        ]
        frame["types"] = [
            get_random_name() for _ in range(n_atoms)
        ]
        frame["res_names"] = [
            get_random_name() for _ in range(n_atoms)
        ]
        frame["res_ids"] = np.random.randint(1, high=n_atoms, size=n_atoms)
        frame["charges"] = np.random.rand(n_atoms) * 2. - 1.
        frame["masses"] = np.random.rand(n_atoms) * 50.
        n_bonds  = np.random.randint(1000, 5000)
        bonds = set()
        for _ in range(n_bonds):
            i = np.random.randint(0, n_atoms)
            j = np.random.randint(0, n_atoms)
            if i == j:
                continue
            # it'd only write it once if there are duplicates
            # + we want this to be the same format as it'll be read back
            smaller = min(i,j)
            larger = max(i,j)
            bonds.add((smaller, larger))
        frame["bonds"] = bonds
        # and write it to file

    print("writing it to file")
    for frame_num, frame in tqdm(enumerate(frames), total=len(frames)):
        w.new_frame(
            frame_num,
            frame_num * 5000,
            frame_num * 0.1,
            frame["n_atoms"],
        )
        w.register_frame_atoms(
            frame["names"],
            frame["res_names"],
            frame["res_ids"],
            frame["types"],
            frame["charges"],
            frame["masses"]
        )
        w.register_frame_bonds(
            frame["bonds"]
        )
        w.write_frame()
    w.finish()


    print("reading from file and verifying")
    r = TopTrajReader("out.toptraj")
    assert r.title == b"example title"
    assert r.initial_molecules[0] == (b"a", 1)
    assert r.initial_molecules[1] == (b"b", 2)
    assert r.initial_molecules[2] == (b"c", 3)

    for i, frame in tqdm(enumerate(frames), total=len(frames)):
        f = r.read_frame()
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
        assert len(f.bonds) == len(frame["bonds"])
        assert len(set(f.bonds) - frame["bonds"]) == 0
        assert len(frame["bonds"] - set(f.bonds)) == 0

    assert r.read_frame() is None


if __name__ == "__main__":
    test_toptraj_writer()
