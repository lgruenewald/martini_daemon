from martini_daemon import PeriodicBox
import numpy as np


def test_which_atoms_within_distance():
    pbc = PeriodicBox.cubic(5.0)
    n_atoms = 50000
    r = 0.5

    atoms = np.random.rand(n_atoms * 3).reshape((n_atoms, 3)) * 5.0

    slow = set()
    for i, atom in enumerate(atoms):
        if pbc.distance(atom, atoms[0]) < r:
            slow.add(i)

    fast = pbc.which_atoms_within_distance(atoms, {0}, r)

    assert fast == slow


def test_which_atoms_within_distance_multiple():
    pbc = PeriodicBox.cubic(5.0)
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


# TODO rest of the stuff
