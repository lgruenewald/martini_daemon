#!/usr/bin/env python3
"""
topology.py

TLDR: - topology is a list of fragment types and reaction templates
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


@dataclass
class ReactionTemplate:
    name: str
    r1: str
    r2: str
    p1: str
    distance_max: float


@dataclass
class Topology:
    frag_fragments: list[FragFragment]
    mol_fragments: list[MolFragment]
    reaction_list: list[ReactionTemplate]
