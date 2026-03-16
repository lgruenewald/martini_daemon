cpdef double pdist(double[:] v1, double[:] v2, double[:] size)
cpdef double pcos_angle(double[:] p1, double[:] p2, double[:] p3, double[:] box)
cpdef pdihedral(double[:] p4, double[:] p3, double[:] p2, double[:] p1, double[:] box)
cpdef bint cross_box(double[:] v1, double[:] v2, double[:] size)

ctypedef long long i64
