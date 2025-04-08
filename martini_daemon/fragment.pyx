from enum import Enum
from typing import Any
from .utils cimport i64
import numpy as np

cdef class Fragment():
    cdef str name
    cdef i64[:] atoms
    cdef i64 frag_id
    graph: GraphFragment

    def __init__(self, graph, i64 id, i64 n_particles):
        self.name = graph.name
        self.atoms = np.zeros((n_particles,))
        self.frag_id = id
        self.graph = graph

    cpdef i64 index_atom(self, i64 i):
        return self.atoms[i]