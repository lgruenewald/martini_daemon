#!/usr/bin/env python

import os
import numpy as np
from math import isclose, pi, cos, tau


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


def pdihedral(p1, p2, p3, p4, box):
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


def test_utils():
    box = np.array([5., 5., 5.])
    box2 = np.array([1., 2., 3.])

    origin = np.array([0., 0., 0.])
    v1 = np.array([0., 0., 1.])
    v2 = np.array([0., 0., 4.])
    v3 = np.array([1., 0., 0.])
    v4 = np.array([2.5, 0., 0.])
    v5 = np.array([1., 0., 1.])

    # testing pdist
    assert isclose(pdist(v1, v2, box), 2.)
    assert isclose(pdist(v1, v3, box), np.sqrt(2))
    assert isclose(pdist(v3, v4, box2), 0.5)

    # testing crossbox
    assert cross_box(v1, v2, box)
    assert not cross_box(v1, v3, box)

    # testing pcos_angle
    assert isclose(pcos_angle(v1, origin, v3, box), 0.)  # 90 degrees
    assert isclose(pcos_angle(v1, origin, v2, box), -1.)  # 180 degrees
    assert isclose(pcos_angle(v1, origin, v5, box), cos(pi / 4))  # 45 degrees

    # corners of a tetrahedron constructed from the box centered at origin,
    # edge length 1, but a random multiple of 5 added to coordinates
    # t1 to t4 are essentially
    # 1/2 1/2 1/2
    # -1/2 -1/2 1/2
    # -1/2 1/2 -1/2
    # 1/2 -1/2 -1/2
    t1 = np.array([5.5, 0.5, 0.5])
    t2 = np.array([-10.5, 9.5, -4.5])
    t3 = np.array([4.5, 5.5, -0.5])
    t4 = np.array([0.5, -0.5, -0.5])
    assert isclose(pcos_angle(t1, t2, t3, box), cos(pi / 3))  # 60 degrees
    assert isclose(pcos_angle(t2, t3, t4, box), cos(pi / 3))  # 60 degrees
    assert isclose(pcos_angle(t1, t4, t3, box), cos(pi / 3))  # 60 degrees
    assert isclose(pcos_angle(t1, t2, t4, box), cos(pi / 3))  # 60 degrees

    # does a all 4 points on a single line type of dihedral work?
    l1 = np.array([1., 0., 0.])
    l2 = np.array([2., 0., 0.])
    l3 = np.array([3., 0., 0.])
    l4 = np.array([4., 0., 0.])
    assert isclose(pdihedral(l1, l2, l3, l4, box), 0.)  # 0 degrees
    assert isclose(pdihedral(l1, l2, l4, l3, box), 0.)  # 0 degrees
    
    # testing dihedrals
    # first 3 have random multiples of the box size added to them
    # but are 1 0 0
    #         2 0 1
    #         3 0 1
    # essentially
    d1 = np.array([1., 5., 0.])
    d2 = np.array([7., 0., 6.])
    d3 = np.array([8., 5., 1.])
    # the dihedral angle 0
    d4_0 = np.array([4., 0., 0.])
    # dihedral angle 45 degrees
    d4_45 = np.array([4., 1.0, 0.])
    # 90 degrees
    d4_90 = np.array([4., 1.0, 1.0])
    # 135 degrees
    d4_135 = np.array([4., 1., 2.])
    # 180
    d4_180 = np.array([4., 0., 2.])
    # 180 + 45 = 225 degrees
    d4_225 = np.array([4., -1., 2.])
    # 270 degrees
    d4_270 = np.array([4., -1., 1.])

    assert isclose(pdihedral(d1, d2, d3, d4_0, box), 0)
    assert isclose(pdihedral(d1, d2, d3, d4_45, box), pi / 4)
    assert isclose(pdihedral(d1, d2, d3, d4_90, box), pi / 2)
    assert isclose(pdihedral(d1, d2, d3, d4_135, box), 3 * pi / 4)
    assert isclose(pdihedral(d1, d2, d3, d4_180, box), pi)
    assert isclose(pdihedral(d1, d2, d3, d4_225, box), 5 * pi / 4)
    assert isclose(pdihedral(d1, d2, d3, d4_270, box), 3 * pi / 2)

    return True


if __name__ == "__main__":
    test_utils()
