# detection algo, attempt 2
import numpy as np
from .reaction_template import ReactionTemplate


def pbc_diff(p1, p2, size):
    """
    Difference (p1-p2) of two positions in a pbc box of size.
    """
    diff = p1 - p2
    base = np.floor(diff / size + 0.5) * size
    return diff - base


def pbc_dist(p1, p2, size):
    """
    Distance between two points in a pbc box of size.
    """
    return np.sqrt(np.sum(pbc_diff(p1, p2, size)**2, axis=-1))


def pbc_angle_cos(p1, p2, p3, size):
    """
    Cosine of the angle enclosed by p1, p2 and p3, p2 being the central atom,
    in a pbc box of size.
    """
    v1 = pbc_diff(p2, p1, size)
    v2 = pbc_diff(p2, p3, size)
    return np.dot(v1, v2) / np.sqrt(v1.dot(v1)) / np.sqrt(v2.dot(v2))
