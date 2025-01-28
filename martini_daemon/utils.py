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


"""
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double pdist(double[:] v1, double[:] v2, double[:] size):
    # periodic boundary adjusted difference between two Vec3
    # based on ReferenceForce.cpp getDeltaRPeriodic from OpenMM

    cdef double diff
    cdef double base
    cdef double sum = 0
    for i in range(3):
        diff = v1[i] - v2[i]
        base = floor(diff / size[i] + 0.5) * size[i]
        sum += pow(diff-base, 2)
    return sqrt(sum)
"""


def pdist(v1, v2, size):
    diff = v1 - v2
    base = np.floor(diff / size + 0.5) * size
    return np.sqrt(np.sum((diff - base)**2))


def cross_box(v1, v2, size):
    dist = np.sqrt(np.sum((v1 - v2)**2))
    pbc_dist = pdist(v1, v2, size)
    if np.abs(dist - pbc_dist) > 0.00001:
        return True


def test_utils():
    example_box = np.array([5., 5., 5.])
    example_box_2 = np.array([1., 2., 3.])

    assert (np.abs(pdist(np.array([0., 0., 1.]), np.array([0., 0., 4.]), example_box) - 2.) < 0.00001)
    assert (np.abs(pdist(np.array([0., 0., 1.]), np.array([1., 0., 0.]), example_box) - np.sqrt(2)) < 0.00001)

    assert (np.abs(pdist(np.array([1., 0., 0.]), np.array([2.5, 0., 0.]), example_box_2) - 0.5) < 0.00001)

    return True
