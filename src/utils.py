#!/usr/bin/env python3

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


def pdist(v1, v2, size):
    # periodic boundary adjusted difference between two Vec3
    # based on ReferenceForce.cpp getDeltaRPeriodic from OpenMM
    # size can be a float or a numpy array of length 3
    # TODO great candidate for optimization

    diff = v1 - v2
    base = np.floor(diff / size + 0.5) * size
    return np.sqrt(np.dot(diff-base, diff-base))


example_box = np.array([5., 5., 5.])
example_box_2 = np.array([1., 2., 3.])

assert (np.abs(pdist(np.array([0., 0., 1.]), np.array([0., 0., 4.]), 5.) -
               2.) < 0.00001)
assert (np.abs(pdist(np.array([0., 0., 1.]), np.array([1., 0., 0.]), example_box) -
               np.sqrt(2)) < 0.00001)

assert (np.abs(pdist(np.array([1., 0., 0.]), np.array([2.5, 0., 0.]), example_box_2) - 0.5) < 0.00001)
