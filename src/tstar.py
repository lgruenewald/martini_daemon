#!/usr/bin/env python3
"""
tstar.py

TLDR: - topology is a list of fragment types, reaction templates and fragments
      - fragments are a group of atoms and their internal interactions
      - reaction templates are fragment based rules used by the D/M algorithm
      to do reactions

topology is an object that contains the following information:
- list of fragment types
    - can be molecules (with bond, angle, etc. strengths), all the information
        needed to instantiate a new molecule in the system basically, except
        the coordinates of particles (similar to a single molecule in itp)
    - can be frag's, which are just mappings to define "subfragments" within
        a molecule
- list of reaction templates
    - reaction name, reactants, products, distance cutoff, constraints,
        initator atom indicies
    - no longer contains a "remapping" of atoms from product->reactant
        as that can be achieved by just using a .frag fragment to "remap" and
        then just use the .frag fragment as the reactant rather than .itp one
- list of fragments
    - list of fragments currently in the system
    - "defrag list" - backmapping of particles to fragments that contain them
"""
from dataclasses import dataclass


@dataclass
class FragFragment:
    mol: str
    mame: str
    # list[(in_frag_id, in_itp_id, type, name, is_edge)]
    atoms: list[(int, int, str, str, bool)]


@dataclass
class MolFragment:
    molecule_name: str
    complete: bool  # if true it has no edge atoms
    # atoms: type, resnum, resname, atomname, chargegr, charge, mass
    atoms: list[(str, int, str, str, int, float, float)]
    # harmonic_bonds: i j length force
    harmonic_bonds: list[(int, int, float, float)]
    # harmonic_angles: i j k theta force
    harmonic_angles: list[(int, int, int, float, float)]
    # proper_dihedrals: i j k l theta force mult
    proper_dihedrals: list
    # improper_dihedrals: i j k l theta force
    improper_dihedrals: list
    # exclusions: list of id's
    exclusions: list[list[int]]
    # constraints: i j length
    constraints: list[(int, int, float)]

    def add_exclusion(self, i, j):
        """Adds an exclusion to the list of exclusions
        Does not add duplicate exclusions.
        Adds it for both the i->j and j->i directions
        Bounds checks and adds empty lists as needed"""

        bigger = max(i, j)
        if bigger >= len(self.atoms):
            raise ValueError(f"{bigger} is out of bounds. "
                             "Particle count in MolFragment: {len(self.atoms)}"
                             ", while add_exclusion({i}, {j}) was called.")
        while len(self.exclusions) - 1 > bigger:
            self.exclusions.append([])

        # this is a slow algorithm, but simple
        if i not in self.exclusions[j]:
            self.exclusions[j].append(i)
            if j in self.exclusions[i]:
                # should never happen actually if this algorithm behaves as
                # i expect it to
                raise AssertionError("ij desynced")
            self.exclusions[i].append(j)

        if j not in self.exclusions[i]:
            raise AssertionError("ij desynced")


@dataclass
class ReactionTemplate:
    name: str
    r1: str
    r2: str
    p1: str
    distance_max: float


class Fragment():
    """Helper class for building and storing fragment information
    fragment name is used for its type
    fragment name can be molname (from .itps) or fragname (from .frag)

    Constructed by either DaemonTopology (when instantiating fragments)

    or constructed by the D/M algorithm (when dynamically generating complete
    fragments)
    """

    name: str
    bonds: list[int]
    angles: list[int]
    dihedrals: list[int]
    impropers: list[int]
    exclusions: list[int]
    constraints: list[int]

    def __init__(self, name):
        self.name = name
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.impropers = []
        self.exclusions = []
        self.constraints = []


class TopStar():
    # T* fragment and defrag list
    frag_list: list[Fragment]
    # for every part_id have a list of fragments it is in
    # type: list[list[(frag_id, in_fragment_id)]]
    # this list has to be at least 1 long per particle, and the first element
    # should always be the reference to the complete fragment
    defrag_list: list[list[(int, int)]]

    # Formerly DaemonTopology:
    frag_fragments: list[FragFragment]
    mol_fragments: list[MolFragment]
    reaction_list: list[ReactionTemplate]
