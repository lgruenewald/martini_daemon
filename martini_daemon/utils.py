import os
import numpy as np


def backup_try(path):
    if os.path.isfile(path):
        bkup_num = 2
        bkup_path = f"#{path}.1#"
        while os.path.isfile(bkup_path):
            bkup_num += 1
            bkup_path = f"#{path}.{bkup_num}#"
        os.rename(path, bkup_path)
        print(f"Backed up {path} to {bkup_path}")


def psub(v1, v2, size):
    diff = v1 - v2
    base = np.floor(diff / size + 0.5) * size
    return diff - base


def pdist(v1, v2, size):
    pdiff = psub(v1, v2, size)
    return np.sqrt(np.sum(pdiff**2))


def cross_box(v1, v2, size):
    return np.any((v1 - v2) / size > 0.5)
#    dist = np.sqrt(np.sum((v1 - v2)**2))
#    pbc_dist = pdist(v1, v2, size)
#    if np.abs(dist - pbc_dist) > 0.00001:
#        return True


def test_utils():
    example_box = np.array([5., 5., 5.])
    example_box_2 = np.array([1., 2., 3.])

    assert (np.abs(pdist(np.array([0., 0., 1.]), np.array([0., 0., 4.]), example_box) - 2.) < 0.00001)
    assert (np.abs(pdist(np.array([0., 0., 1.]), np.array([1., 0., 0.]), example_box) - np.sqrt(2)) < 0.00001)

    assert (np.abs(pdist(np.array([1., 0., 0.]), np.array([2.5, 0., 0.]), example_box_2) - 0.5) < 0.00001)

    return True
