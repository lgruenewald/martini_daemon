#!/usr/bin/env python3

class Fragment():
    """Helper class for building and storing fragment information
    fragment name is used for its type
    fragment name can be molname (from .itps) or fragname (from .frag)
    """

    def __init__(self, topstar, name):
        self.name = name
        self._topstar = topstar
        self._bonds = []
        # list of (part_id, is_edge:bool)
        self._particles = []

    def add_particle(self, ptype, is_edge, frag_id):
        # FIXME frag_id can be moved probably
        i = self._topstar.part_list.add(ptype)
        self.add_existing_particle(i, is_edge, frag_id)

    def add_existing_particle(self, id, is_edge, frag_id):
        # FIXME frag_id can be moved probably
        self._particles.append((id, is_edge))
        in_frag_id = len(self._particles) - 1
        self._topstar.part_list.add_defrag((frag_id, in_frag_id))

    def add_bond(self, mid1, mid2, funct, k, l0):
        # called with in#include_fragment_id
        # (so not the indices in topstar.partlist)
        if mid1 >= len(self._particles) or mid2 >= len(self._particles):
            raise ValueError("add_bond should be called with per molecule"
                             "particle indices")
        self._bonds.append((mid1, mid2, funct, k, l0))
        id1 = self._particles[mid1]
        id2 = self._particles[mid2]
        bond_id = self._topstar.bondlist.add(id1, id2, funct, k, l0)
        self._bonds.append(bond_id)

    def add_existing_bond(self, id):
        self._bonds.append(id)


class TopStar():
    frag_list: list[Fragment]
    # for every part_id have a list of fragments it is in
    # type: list[list[(frag_id, in_fragment_id)]]
    # this list has to be at least 1 long per particle, and the first element
    # should always be the reference to the complete fragment
    defrag_list: list[list[(int, int)]]
