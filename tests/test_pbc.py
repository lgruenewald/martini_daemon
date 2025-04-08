import numpy as np
from math import cos, pi, isclose
from martini_daemon.utils import cross_box, pcos_angle, pdihedral, pdist


def test_pbc():
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

    # testing crossbox
    assert cross_box(v1, v2, box)
    assert not cross_box(v1, v3, box)
