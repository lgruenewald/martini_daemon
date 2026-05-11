"""Test the class PeriodicBox."""

import mdtraj
import numpy as np
import pytest

from martini_daemon import PeriodicBox

boxes = [
    PeriodicBox.cubic(10.0),
    PeriodicBox.orthogonal(5.0, 7.0, 9.0),
    PeriodicBox.orthogonal(2.0, 5.0, 9.0),
    PeriodicBox([10.0, 0.0, 0.0], [2.4, 10.0, 0.0], [2.4, -4.0, 10.0]),
    PeriodicBox([10.0, 0.0, 0.0], [-4.0, 10.0, 0.0], [-4.0, -4.0, 10.0]),
    PeriodicBox([10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [2.0, 4.0, 10.0]),
    PeriodicBox([10.0, 0.0, 0.0], [0.0, 5.0, 0.0], [2.0, 2.0, 2.0]),
]


@pytest.mark.parametrize("pbc", boxes)
def test_distances(pbc: PeriodicBox) -> None:
    """Test PeriodicBox.distance() and PeriodicBox.distance_squared()."""
    n_atoms = 100
    xyz = np.random.rand(n_atoms * 3).reshape((1, n_atoms, 3)) * 100.0 - 50.0
    traj = mdtraj.Trajectory(xyz, None)
    traj.unitcell_vectors = np.array([[pbc.a, pbc.b, pbc.c]])
    for i in range(n_atoms):
        for j in range(i):
            dist_ref = float(mdtraj.compute_distances(traj, [(i, j)])[0][0])
            dist = pbc.distance(xyz[0, i], xyz[0, j])
            dist_sq = pbc.distance_squared(xyz[0, i], xyz[0, j])
            # 1e-3 because mdtraj is f32, and that's the precision of .gro/.xtc files anyway, so should be enough
            assert np.isclose(dist_ref, dist, atol=1e-3), (
                f"Bad distance. Coords: {xyz[0, i]}, {xyz[0, j]}.\n"
                + f"Box: {pbc.a} {pbc.b} {pbc.c}.\n"
                + f"Distance martini_daemon: {dist} nm.\n"
                + f"Distance mdtraj: {dist_ref} nm.\n"
            )
            assert np.isclose(dist_ref**2, dist_sq, atol=1e-3)


@pytest.mark.parametrize("pbc", boxes)
def test_move_to(pbc: PeriodicBox) -> None:
    """Test PeriodicBox.move_to()."""
    n_atoms = 100
    xyz = np.random.rand(n_atoms * 3).reshape((1, n_atoms, 3)) * 100.0 - 50.0
    traj = mdtraj.Trajectory(xyz, None)
    traj.unitcell_vectors = np.array([[pbc.a, pbc.b, pbc.c]])
    for i in range(n_atoms):
        for j in range(i):
            dist_ref = float(mdtraj.compute_distances(traj, [(i, j)])[0][0])
            moved_j = np.array(pbc.move_to(xyz[0, i], xyz[0, j]))
            dist = np.linalg.norm(xyz[0, i] - moved_j)
            assert np.isclose(dist_ref, dist, atol=1e-3), (
                f"Bad move_to. Ref I: {xyz[0, i]} J: {xyz[0, j]}.\n"
                f"Moved J: {moved_j}.\n"
                + f"Box: {pbc.a} {pbc.b} {pbc.c}.\n"
                + f"Distance martini_daemon (move_to): {dist} nm.\n"
                + f"Distance mdtraj: {dist_ref} nm.\n"
            )


@pytest.mark.parametrize("pbc", boxes)
def test_angles(pbc: PeriodicBox) -> None:
    """Test PeriodicBox.angle() and PeriodicBox.cos_angle()."""
    n_atoms = 50
    xyz = np.random.rand(n_atoms * 3).reshape((1, n_atoms, 3)) * 100.0 - 50.0
    traj = mdtraj.Trajectory(xyz, None)
    traj.unitcell_vectors = np.array([[pbc.a, pbc.b, pbc.c]])
    for i in range(n_atoms):
        for j in range(i):
            for k in range(j):
                angle_ref = float(mdtraj.compute_angles(traj, [(i, j, k)])[0][0])
                angle = pbc.angle(xyz[0, i], xyz[0, j], xyz[0, k])
                cos_angle = pbc.cos_angle(xyz[0, i], xyz[0, j], xyz[0, k])
                assert np.isclose(angle_ref, angle, atol=1e-3), (
                    f"Bad angle. pos: {xyz[0, i]} {xyz[0, j]} {xyz[0, k]}\n"
                    + f"daemon: {angle}\n"
                    + f"mdtraj: {angle_ref}\n"
                )
                # tol of ~0.01 deg
                assert np.isclose(np.cos(angle_ref), cos_angle, atol=1e-3)


@pytest.mark.parametrize("pbc", boxes)
def test_dihedrals(pbc: PeriodicBox) -> None:
    """Test PeriodicBox.dihedral()."""
    n_atoms = 10
    xyz = np.random.rand(n_atoms * 3).reshape((1, n_atoms, 3)) * 100.0 - 50.0
    traj = mdtraj.Trajectory(xyz, None)
    # mdtraj raises a performance warning if not f32
    traj.unitcell_vectors = np.array([[pbc.a, pbc.b, pbc.c]], dtype=np.float32)
    for i in range(n_atoms):
        for j in range(i):
            for k in range(j):
                for m in range(k):
                    dih_ref = float(
                        mdtraj.compute_dihedrals(traj, [(i, j, k, m)])[0][0]
                    )
                    dih = pbc.dihedral(xyz[0, i], xyz[0, j], xyz[0, k], xyz[0, m])
                    # tol of ~0.1 deg
                    # unfortunately rarely errors of ~0.02 deg seem to occur
                    assert np.isclose(dih, dih_ref, atol=1e-2)


@pytest.mark.parametrize("pbc", boxes)
def test_move_within(pbc: PeriodicBox) -> None:
    """Test PeriodicBox.move_within()."""
    n_atoms = 100
    xyz = np.random.rand(n_atoms * 3).reshape((1, n_atoms, 3)) * 100.0 - 50.0
    traj = mdtraj.Trajectory(xyz, None)
    traj.unitcell_vectors = np.array([[pbc.a, pbc.b, pbc.c]])
    for i in range(n_atoms):
        within = pbc.move_within(xyz[0, i])
        assert 0 < within[0] < pbc.a[0]
        assert 0 < within[1] < pbc.b[1]
        assert 0 < within[2] < pbc.c[2]


@pytest.mark.parametrize("pbc", boxes)
def test_which_atoms_within_distance(pbc: PeriodicBox) -> None:
    """Test PeriodicBox.which_atoms_within_distance()."""
    n_atoms = 50000
    r = 0.5

    atoms = np.random.rand(n_atoms * 3).reshape((n_atoms, 3)) * 10.0

    slow = set()
    for i, atom in enumerate(atoms):
        if pbc.distance(atom, atoms[0]) < r:
            slow.add(i)

    fast = pbc.which_atoms_within_distance(atoms, {0}, r)

    assert fast == slow, (
        "Test which atoms within distance fail.\n"
        + f"Set slow size: {len(slow)}, set fast size: {len(fast)}\n"
        + f"Diff slow-fast: {slow - fast}\n"
        + f"Diff fast-slow: {fast - slow}\n"
    )


@pytest.mark.parametrize("pbc", boxes)
def test_which_atoms_within_distance_multiple(pbc: PeriodicBox) -> None:
    """Test PeriodicBox.which_atoms_within_distance() with multiple atoms."""
    n_atoms = 500000
    r = 0.5

    atoms = np.random.rand(n_atoms * 3).reshape((n_atoms, 3)) * 5.0

    selection = set(np.random.randint(n_atoms, size=10))

    slow = set()
    for i, atom in enumerate(atoms):
        for s in selection:
            if pbc.distance(atom, atoms[s]) < r:
                slow.add(i)
                break

    fast = pbc.which_atoms_within_distance(atoms, selection, r)

    assert fast == slow
