#!/usr/bin/env python3
"""
topstar.py

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
from sysstar import SysStar

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

        if i == j:
            raise ValueError(f"add_exclusion({i}, {j}) where i==j was called.")
        bigger = max(i, j)
        if bigger >= len(self.atoms):
            raise ValueError(f"{bigger} is out of bounds. "
                             f"Particle count in MolFragment: {len(self.atoms)}"
                             f", while add_exclusion({i}, {j}) was called.")
        while len(self.exclusions) <= bigger:
            self.exclusions.append([])

        # this is a slow algorithm, but simple
        if i not in self.exclusions[j]:
            self.exclusions[j].append(i)
            if j in self.exclusions[i]:
                # should never happen actually if this algorithm behaves as
                # i expect it to
                raise AssertionError(f"ij desynced, {j} "
                                     f"already in self.exclusions[{i}]")
            self.exclusions[i].append(j)

        if j not in self.exclusions[i]:
            raise AssertionError(f"ij desynced, {j} not in "
                                 f"self.exclusions[{i}]")


@dataclass
class ReactionTemplate:
    name: str
    r1: str
    r2: str
    p1: str
    distance_max: float


@dataclass
class Fragment():
    """Helper class for building and storing fragment information
    fragment name is used for its type
    fragment name can be molname (from .itps) or fragname (from .frag)

    Constructed by either DaemonTopology (when instantiating fragments)

    or constructed by the D/M algorithm (when dynamically generating complete
    fragments)
    """

    name: str
    particles: list[int]
    bonds: list[int]
    angles: list[int]
    dihedrals: list[int]
    impropers: list[int]
    exclusions: list[int]
    constraints: list[int]

    def __init__(self, name):
        self.name = name
        self.particles = []
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
    type_lookup: dict[str, FragFragment | MolFragment]
    reaction_list: list[ReactionTemplate]

    system: SysStar

    def __init__(self, system):
        self.frag_list = []
        self.defrag_list = []
        self.frag_fragments = []
        self.mol_fragments = []
        self.reaction_list = []
        self.type_lookup = {}
        self.system = system

    def new_mol_fragment(self, name: str):
        if self.type_lookup.get(name):
            raise ValueError(f"Second definition of fragment type {name}")
        mol_fragment = MolFragment(name, True, [], [], [], [], [], [], [])
        self.mol_fragments.append(mol_fragment)
        self.type_lookup[name] = mol_fragment

    def instantiate(self, frag_name):
        """Takes a name of a mol fragment, creates new particles for it in
        the system and the corresponding interactions as well.

        Later, when developing the D/M algorithm a modified version of this
        should be created for reusing existing particles.
        """

        frag = self.type_lookup.get(frag_name)  # frag type
        if frag is None:
            raise ValueError(f"Can't find mol {frag_name}")
        elif not isinstance(frag, MolFragment):
            raise ValueError(f"Attempt to instantiate {frag_name}, but it's"
                             " not a mol fragment type."
                             f" It is: {frag}")
        inst = Fragment(frag_name)
        self.frag_list.append(inst)
        # particles
        index0 = self.system.len_particles()
        for atom in frag.atoms:
            type, resnum, resname, atomname, chargegr, charge, mass = \
                atom
            # TODO if mass is -1 use a default
            p = self.system.add_particle(atomname, type, charge, mass)
            inst.particles.append(p)
        # bonds
        for bond in frag.harmonic_bonds:
            i, j, length, force = bond
            b = self.system.add_bond(index0 + i, index0 + j, length, force)
            inst.bonds.append(b)
        # angles
        for angle in frag.harmonic_angles:
            i, j, k, theta, force = angle
            a = self.system.add_angle(i + index0, j + index0, k + index0,
                                       theta, force)
            inst.angles.append(a)
        # proper dihedrals
        for dih in frag.proper_dihedrals:
            i, j, k, l, theta, force, mult = dih
            d = self.system.add_proper_dihedral(i + index0, j + index0, k + index0,
                                          l + index0, theta, force, mult)
            inst.dihedrals.append(d)
        # improper dihedrals
        for imp in frag.improper_dihedrals:
            i, j, k, l, theta, force = imp
            d = self.system.add_improper_dihedral(i + index0, j + index0, k + index0,
                                          l + index0, theta, force)
            inst.impropers.append(d)
        # exclusions
        for i, excl in enumerate(frag.exclusions):
            for j in excl:
                if i < j:
                    e = self.system.add_exclusion(i + index0, j + index0)
                    inst.exclusions.append(e)
        # constraints
        for cons in frag.constraints:
            i, j, length = cons
            c = self.system.add_constraint(i + index0, j + index0, length)
            inst.constraints.append(c)
