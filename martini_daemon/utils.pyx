import os
import numpy as np
from math import tau
cimport cython
cimport numpy as cnp
cdef extern from "math.h":
    double floor(double v)
    double sqrt(double v)
    double atan2(double u, double v)
    double fabs(double v)

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


def format_time(total):
    c = int(total)
    seconds = c % 60
    c //= 60
    minutes = c % 60
    c //= 60
    hours = c % 24
    c //= 24
    days = c

    res = f"{hours:02}:{minutes:02}:{seconds:02}"
    if days > 0:
        res = f"{days}d " + res
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
cdef cnp.ndarray[cnp.float64_t] psub(double[:] v1, double[:] v2, double[:] size):
    cdef cnp.ndarray[cnp.float64_t] res = np.zeros((3,), dtype=np.float64)
    cdef double diff
    cdef double base
    for i in range(3):
        diff = v1[i] - v2[i]
        base = floor(diff / size[i] + 0.5) * size[i]
        res[i] = diff - base
    return res


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double pdist(double[:] v1, double[:] v2, double[:] size):
    cdef double[:] pdiff = psub(v1, v2, size)
    cdef double sum = 0.0
    for i in range(3):
        sum += pdiff[i] * pdiff[i]
    return sqrt(sum)

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double pcos_angle(double[:] p1, double[:] p2, double[:] p3, double[:] box):
    cdef double[:] v1 = psub(p2, p1, box)
    cdef double[:] v2 = psub(p2, p3, box)
    cdef double v1_dot_v2 = 0.0
    cdef double v1_dot_v1 = 0.0
    cdef double v2_dot_v2 = 0.0
    for i in range(3):
        v1_dot_v1 += v1[i] * v1[i]
        v1_dot_v2 += v1[i] * v2[i]
        v2_dot_v2 += v2[i] * v2[i]
    return v1_dot_v2 / sqrt(v1_dot_v1) / sqrt(v2_dot_v2)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef cnp.ndarray[cnp.float64_t] cross_product(u: double[:], v: double[:]):
    cdef cnp.ndarray[cnp.float64_t] res = np.zeros((3,), dtype=np.float64)
    res[0] = u[1] * v[2] - u[2] * v[1]
    res[1] = u[2] * v[0] - u[0] * v[2]
    res[2] = u[0] * v[1] - u[1] * v[0]
    return res


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef pdihedral(double[:] p4, double[:] p3, double[:] p2, double[:] p1, double[:] box):
    # difference vectors
    cdef double[:] u1 = psub(p2, p1, box)
    cdef double[:] u2 = psub(p3, p2, box)
    cdef double[:] u3 = psub(p4, p3, box)
    # cross products
    cdef double[:] v1 = cross_product(u1, u2)
    cdef double[:] v2 = cross_product(u2, u3)
    cdef double[:] cross_v1_v2 = cross_product(v1, v2)
    # x, y coords (cos and sin of the angle)
    cdef double x = 0.0
    cdef double y = 0.0
    cdef double norm_u2 = sqrt(u2[0] * u2[0] + u2[1] * u2[1] + u2[2] * u2[2])
    for i in range(3):
        x += v1[i] * v2[i]  # v1 . v2
        y += u2[i] * cross_v1_v2[i] / norm_u2  # u2 . (v1 x v2) / |u2|
    cdef double res = atan2(y, x)
    if res < 0.:
        return res + tau
    else:
        return res

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef bint cross_box(double[:] v1, double[:] v2, double[:] size):
    for i in range(3):
        if fabs((v1[i] - v2[i]) / size[i]) > 0.5:
            return True
    return False
