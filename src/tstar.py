#!/usr/bin/env python3

import typing

TopStar = typing.NewType("TopStar", None)


# ==== System layer ====
class ParticleList():
    """All particles in T*
    indices are called part_id
    values are (part_name: string, part_type: string, list(frag_id, in_frag_id))
    """
    _topstar: TopStar
    _part_list: list[(int, list[int, int])]

    def __init__(self, topstar):
        self._part_list = []
        self._topstar = topstar

    def add(self, part_name, part_type):
        """Adds a particle to the list, and returns its part_id"""
        self._part_list.append((part_name, part_type, []))
        return len(self._part_list) - 1

    def add_defrag(self, part_id, frag_id, in_frag_id):
        if part_id >= len(self._part_list):
            raise ValueError("Invalid index")
        self._part_list[part_id][1].append((frag_id, in_frag_id))


class BondList():
    """All bonds in T*
    indices are called bond_id
    values are (i: part_id, j: part_id, funct: int, k: float, length: float)
    """
    _topstar: TopStar
    _bond_list: list[(int, int, int, float, float)]

    def __init__(self, topstar):
        self._bond_list = []
        self._topstar = topstar

    def add(self, part_id_i, part_id_j, funct, k, length):
        """Adds a bond to the list, returns its bond_id"""
        self._bond_list.append((part_id_i, part_id_j, funct, k, length))
        return len(self._bond_list) - 1


# ==== Fragment layer ====
class Fragment():
    """Helper class for building and storing fragment information
    fragment name is used for its type
    fragment name can be molname (from .itps) or fragname (from .frag)
    """
    _topstar: TopStar

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
        # called with in_fragment_id
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


class FragmentList():
    _topstar: TopStar
    _frag_list: list[Fragment]
    # for every part_id have a list of fragments it is in
    # type: list[list[(frag_id, in_fragment_id)]]
    _defrag_list: list[list[(int, int)]]


# ==== Top Star ====
class TopStar():
    # S*
    part_list: ParticleList
    bond_list: BondList
    # T*
    frag_list: FragmentList

    def __init__(self):
        self.part_list = ParticleList()
        self.bond_list = BondList()
        self.frag_list = FragmentList()


if __name__ == "__main__":
    a = ParticleList()
    print(a.list)
    a.add_particle(13)
    print(a.list)
