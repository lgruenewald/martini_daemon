#!/usr/bin/env python

import os
import numpy as np
from math import tau


def backup_try(path):
    parent, filename = os.path.split(path)
    if os.path.isfile(path):
        bkup_num = 0
        bkup_path = path
        while os.path.isfile(bkup_path):
            bkup_num += 1
            bkup_path = os.path.join(parent, f"#{filename}.{bkup_num}#")
        os.rename(path, bkup_path)
        print(f"Backed up {path} to {bkup_path}")


def psub(v1, v2, size):
    diff = v1 - v2
    base = np.floor(diff / size + 0.5) * size
    return diff - base


def pdist(v1, v2, size):
    pdiff = psub(v1, v2, size)
    return np.sqrt(np.sum(pdiff**2))


def pcos_angle(p1, p2, p3, box):
    """
        Returns the cosine of the angle, with p2 as the central atom
        in a periodic boundary system defined by box
    """
    # the two vectors going from central atom (2) to the edge atoms
    # aware of periodic boxes
    v1 = psub(p2, p1, box)
    v2 = psub(p2, p3, box)
    # to get the cosine of the angle divide the dot product by the
    # lengths of the vectors
    return np.dot(v1, v2) / np.sqrt(v1.dot(v1)) / np.sqrt(v2.dot(v2))


def pdihedral(p4, p3, p2, p1, box):
    """
        Returns the dihedral between the four particles, in radians
    """
    # the three difference vectors
    u1 = psub(p2, p1, box)
    u2 = psub(p3, p2, box)
    u3 = psub(p4, p3, box)
    # cross products
    v1 = np.cross(u1, u2)
    v2 = np.cross(u2, u3)
    # x,y coords (cos and sin of the angle)
    x = np.dot(v1, v2)
    y = np.dot(u2, np.cross(v1, v2)) / np.linalg.norm(u2)
    # cartesian -> polar, but we only need the angle
    res = np.arctan2(y, x)
    return res + tau if res < 0. else res


def cross_box(v1, v2, size):
    return np.any(np.abs((v1 - v2) / size) > 0.5)
