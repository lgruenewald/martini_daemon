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
from fnmatch import fnmatch
import sys


@dataclass
class FragFragment:
    name: str
    mol: str
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
                             f"Particle count in MolFragment:{len(self.atoms)}"
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
    edge: list[bool]
    frag_id: int

    def __init__(self, name, id):
        self.name = name
        self.particles = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.impropers = []
        self.exclusions = []
        self.constraints = []
        self.edge = []
        self.frag_id = id


class TopStar():
    # === Live molecule/fragment information ===
    # T* fragment and defrag list
    frag_list: list[Fragment]
    # for every part_id have a list of fragments it is in
    # type: list[list[(frag_id, in_fragment_id, is_edge)]]
    # this list has to be at least 1 long per particle, and the first element
    # should always be the reference to the complete fragment
    defrag_list: list[list[(int, int, bool)]]

    # === Fragment Types ===
    # type name -> list of subfrag names
    subfrag_map: dict[str, list[str]]

    # type name -> type
    type_lookup: dict[str, FragFragment | MolFragment]

    # === Reaction templates ===
    reaction_list: list[ReactionTemplate]

    reactive_pairs: dict[(str, str), ReactionTemplate]

    reactive_types: set[str]

    system: SysStar

    def __init__(self, system):
        self.frag_list = []
        self.defrag_list = []
        self.frag_fragments = []
        self.mol_fragments = []
        self.reaction_list = []
        self.type_lookup = {}
        self.subfrag_map = {}
        self.reactive_pairs = {}
        self.reactive_types = set()
        self.system = system

    def new_mol_fragment(self, name: str):
        if self.type_lookup.get(name):
            raise ValueError(f"Second definition of fragment type {name}")
        mol_fragment = MolFragment(name, True, [], [], [], [], [], [], [])
        self.type_lookup[name] = mol_fragment

    def new_frag_fragment(self, name: str, parent: str):
        if self.type_lookup.get(name):
            raise ValueError(f"Redefinition of fragment type {name}")
        frag_fragment = FragFragment(name, parent, [])
        self.type_lookup[name] = frag_fragment
        if self.subfrag_map.get(parent):
            self.subfrag_map[parent].append(name)
        else:
            self.subfrag_map[parent] = [name]

    def new_reaction(self, reaction: ReactionTemplate):
        self.reaction_list.append(reaction)
        self.reactive_pairs[(reaction.r1, reaction.r2)] = reaction
        self.reactive_pairs[(reaction.r2, reaction.r1)] = reaction
        self.reactive_types.add(reaction.r1)
        self.reactive_types.add(reaction.r2)

    def instantiate_subfrag(self, subfrag, inst):
        """Instantiates a subfragment subfrag, with forces in S* as seen in
        inst."""
        frag: FragFragment = self.type_lookup.get(subfrag)
        if frag is None:
            raise ValueError(f"Can't find frag {subfrag}")
        elif not isinstance(frag, FragFragment):
            raise ValueError(f"Attempt to recursively instantiate {subfrag}, "
                             "but it's not a frag fragment type."
                             f"It is: {subfrag}")
        subinst_id = len(self.frag_list)
        subinst = Fragment(subfrag, subinst_id)
        self.frag_list.append(subinst)
        normal_indices = set()
        edge_indices = set()
        all_indices = set()
        # particles
        for i, atom in enumerate(frag.atoms):
            # remapping of parent_id's to frag_id's
            frag_id, parent_id, type, name, is_edge = atom
            # frag_id's must be in order
            if frag_id != i:
                raise ValueError("Atom IDs in [ fragatoms ] should be ordered")
            # parent_id must be i range
            if parent_id < 0 or parent_id >= len(inst.particles):
                raise ValueError("Parent particle ID out of range")
            part_index = inst.particles[parent_id]
            # name and type must match S*
            sname, stype = self.system.get_particle_name_type(part_index)
            if not fnmatch(sname, name) or not fnmatch(stype, type):
                raise ValueError("During subfrag instantiation name/type "
                                 f"doesn't match. Expected name {sname} "
                                 f"type {stype}. Got name {name} type {type}.")
            # can only add non edge if it's not an edge in parent
            if not is_edge and inst.edge[parent_id]:
                raise ValueError("Attempt to add non-edge atom in subfrag "
                                 "from an edge atom in parent.")
            # build up these sets that are used for constructing the forces
            # prevent double inclusion of the same atom
            if part_index in all_indices:
                raise ValueError("Double inclusion of atom in subfrag")
            all_indices.add(part_index)
            if is_edge:
                edge_indices.add(part_index)
            else:
                normal_indices.add(part_index)
            # build up subinst
            subinst.particles.append(part_index)
            subinst.edge.append(is_edge)
            # frag_id, in_frag_id, is_edge
            self.defrag_list[part_index].append((subinst_id, i, is_edge))
        # the indices sets must be non overlapping
        if len(normal_indices & edge_indices) != 0:
            print("normal indices: ", normal_indices)
            print("edge: ", edge_indices)
            print("all: ", all_indices)
            raise AssertionError("Normal and edge atoms must not overlap")
        if len(all_indices) != len(normal_indices | edge_indices):
            print("normal indices: ", normal_indices)
            print("edge: ", edge_indices)
            print("all: ", all_indices)
            raise AssertionError("All indices must equal normal | edge")
        # bonds get added if at least one normal atom participates in them
        # bonds to be added must be between atoms inside the subfrag
        # TODO: better error messages for dangling stuff
        # through either analysis of subfrags before instantiation
        # or here through dumping more info
        for bond_id in inst.bonds:
            i, j = self.system.get_bond_members(bond_id)
            if i in normal_indices or j in normal_indices:
                # bond required
                if i not in all_indices or j not in all_indices:
                    raise ValueError("Dangling bond")
                subinst.bonds.append(bond_id)
        # angles get added if the central atom is a normal atom
        # all atoms in such angles must contain only atoms in the subfrag
        for angle_id in inst.angles:
            i, j, k = self.system.get_angle_members(angle_id)
            if j in normal_indices:
                if i not in all_indices or k not in all_indices:
                    raise ValueError("Dangling angle")
                subinst.angles.append(angle_id)
        # dihedrals get added if one of the central atoms is a normal atom
        # all atoms participating must be in the subfrag
        for dih_id in inst.dihedrals:
            i, j, k, l = self.system.get_proper_dihedral_members(dih_id)
            if j not in normal_indices and k not in normal_indices:
                continue
            if i not in all_indices or j not in all_indices or \
                    k not in all_indices or l not in all_indices:
                raise ValueError("Dangling dihedral")
            subinst.dihedrals.append(dih_id)
        # improper dihedrals get added if any of the atoms is a normal atom
        # all atoms participating must be in the subfrag
        for dih_id in inst.impropers:
            i, j, k, l = self.system.get_improper_members(dih_id)
            if i in normal_indices or j in normal_indices or \
                    k in normal_indices or l in normal_indices:
                if i not in all_indices or j not in all_indices or \
                        k not in all_indices or l not in all_indices:
                    raise ValueError("Dangling improper dihedral")
                subinst.impropers.append(dih_id)
        # exclusions, constaints get added if any of the atoms is a normal atom
        # all atoms participating must be in the subfrag
        for excl_id in inst.exclusions:
            i, j = self.system.get_exclusion_members(excl_id)
            if i in normal_indices or j in normal_indices:
                if i not in all_indices or j not in all_indices:
                    raise ValueError("Dangling exclusion")
                subinst.exclusions.append(excl_id)
        for cid in inst.constraints:
            i, j = self.system.get_constraint_members(cid)
            if i in normal_indices or j in normal_indices:
                if i not in all_indices or j not in all_indices:
                    raise ValueError("Dangling constraint")
                subinst.constraints.append(cid)

        # subfrags can contain further subfrags
        self.instantiate_subfrags(subfrag, subinst)

    def instantiate_subfrags(self, frag_name, inst):
        """Instantiates all subfrags of frag_name.
        inst is a Fragment instance that contains the reference to all forces
        in S*."""
        subfrags = self.subfrag_map.get(frag_name)
        if subfrags is not None:
            for subfrag in subfrags:
                self.instantiate_subfrag(subfrag, inst)

    def instantiate(self, frag_name):
        frag = self.type_lookup.get(frag_name)  # frag type
        if frag is None:
            raise ValueError(f"Can't find mol {frag_name}")
        elif not isinstance(frag, MolFragment):
            raise ValueError(f"Attempt to instantiate {frag_name}, but it's"
                             " not a mol fragment type."
                             f" It is: {frag}")
        parts = []
        for in_frag_id, atom in enumerate(frag.atoms):
            type, resnum, resname, atomname, chargegr, charge, mass = atom
            p = self.system.add_particle(atomname, type, charge, mass)
            parts.append(p)
            self.defrag_list.append([])
        return self.instantiate_over_existing(frag_name, parts)

    def instantiate_over_existing(
            self, frag_name, particles, first=False, update=False):
        """Takes a name of a mol fragment, creates new particles for it in
        the system and the corresponding interactions as well.
        Recursively instantiates all subfragments too.
        """

        frag = self.type_lookup.get(frag_name)  # frag type
        if frag is None:
            raise ValueError(f"Can't find mol {frag_name}")
        elif not isinstance(frag, MolFragment):
            raise ValueError(f"Attempt to instantiate {frag_name}, but it's"
                             " not a mol fragment type."
                             f" It is: {frag}")
        inst_id = len(self.frag_list)
        inst = Fragment(frag_name, inst_id)
        self.frag_list.append(inst)
        # particles
        for in_frag_id, part_id in enumerate(particles):
            inst.particles.append(part_id)
            # frag_id, in_frag_id, is_edge
            if first:
                self.defrag_list[part_id].insert(
                    0, (inst_id, in_frag_id, False)
                )
            else:
                self.defrag_list[part_id].append((inst_id, in_frag_id, False))
            inst.edge.append(False)
            if update:
                type, resnum, resname, atomname, chargegr, charge, mass = \
                    frag.atoms[in_frag_id]
                self.system.update_particle(
                    part_id, atomname, type, charge, mass
                )
        # bonds
        for bond in frag.harmonic_bonds:
            i, j, length, force = bond
            b = self.system.add_bond(particles[i], particles[j], length, force)
            inst.bonds.append(b)
        # angles
        for angle in frag.harmonic_angles:
            i, j, k, theta, force = angle
            a = self.system.add_angle(
                particles[i], particles[j], particles[k], theta, force
            )
            inst.angles.append(a)
        # proper dihedrals
        for dih in frag.proper_dihedrals:
            i, j, k, l, theta, force, mult = dih
            d = self.system.add_proper_dihedral(
                particles[i], particles[j], particles[k], particles[l],
                theta, force, mult
            )
            inst.dihedrals.append(d)
        # improper dihedrals
        for imp in frag.improper_dihedrals:
            i, j, k, l, theta, force = imp
            d = self.system.add_improper_dihedral(
                particles[i], particles[j], particles[k], particles[l],
                theta, force
            )
            inst.impropers.append(d)
        # exclusions
        for i, excl in enumerate(frag.exclusions):
            for j in excl:
                if i < j:
                    e = self.system.add_exclusion(particles[i], particles[j])
                    inst.exclusions.append(e)
        # constraints
        for cons in frag.constraints:
            i, j, length = cons
            c = self.system.add_constraint(particles[i], particles[j], length)
            inst.constraints.append(c)

        self.instantiate_subfrags(frag_name, inst)
        return inst

    def destroy_fragment(self, frag: Fragment):
        """args:
        frag: Fragment

        returns: the root molecule fragment from underneath
        """
        # 0. save the molecule fragment
        molfrag_id, _, _ = self.defrag_list[frag.particles[0]][0]
        molfrag = self.frag_list[molfrag_id]

        # 1. remove overlapping fragments
        # get the list of fragments to remove
        for part in frag.particles:
            defrag = self.defrag_list[part]
            i = 0
            while i < len(defrag):
                frag_id, _, is_edge = defrag[i]
                if not is_edge or frag_id == frag.frag_id:
                    # remove frag and frags that contain any of the particles
                    # as normal particles
                    del defrag[i]
                    self.frag_list[frag_id] = None
                else:
                    i += 1

        # 2. remove forces in fragment
        for bond in frag.bonds:
            self.system.remove_bond(bond)
        for angle in frag.angles:
            self.system.remove_angle(angle)
        for dih in frag.dihedrals:
            self.system.remove_proper_dihedral(dih)
        for dih in frag.impropers:
            self.system.remove_improper(dih)
        for excl in frag.exclusions:
            self.system.remove_exclusion(excl)
        for con in frag.constraints:
            # constraints can't be removed, so this will throw an exception
            self.system.remove_constraint(con)

        # 3. return the complete fragment that contained this fragment
        # (can be fragment itself)
        return molfrag

    def remove_fragment(self, frag: Fragment):
        """Removes a fragment from frag_lits and defrag_list
        """
        for part in frag.particles:
            defrag = self.defrag_list[part]
            i = 0
            while i < len(defrag):
                frag_id, _, _ = defrag[i]
                if frag_id == frag.frag_id:
                    del defrag[i]
                else:
                    i += 1
        self.frag_list[frag.frag_id] = None

    def new_dynamic_complete_fragment(self, particles, frag1, frag2, frag_prod,
                                      complete1, complete2):
        # get all forces from frag_prod, complete1, complete2 but
        # exclude those from frag1, frag2 (those were deleted from system)
        # and add a new fragment to the list that contains particles and
        # forces obtained that way
        inst_id = len(self.frag_list)
        inst = Fragment("<dyn>", inst_id)
        self.frag_list.append(inst)
        for i, part in enumerate(particles):
            inst.particles.append(part)
            inst.edge.append(False)
            self.defrag_list[part].insert(0, (inst_id, i, False))

        # multiline editing helps a lot with this, no I didn't type this out
        # by hand FIXME this is getting painful to look at though...
        bonds_yes = set()
        bonds_no = set()
        for bond in frag1.bonds:
            bonds_no.add(bond)
        for bond in frag2.bonds:
            bonds_no.add(bond)
        for bond in frag_prod.bonds:
            bonds_yes.add(bond)
        for bond in complete1.bonds:
            bonds_yes.add(bond)
        for bond in complete2.bonds:
            bonds_yes.add(bond)
        bonds = bonds_yes - bonds_no
        for bond in bonds:
            inst.bonds.append(bond)

        angles_yes = set()
        angles_no = set()
        for angle in frag1.angles:
            angles_no.add(angle)
        for angle in frag2.angles:
            angles_no.add(angle)
        for angle in frag_prod.angles:
            angles_yes.add(angle)
        for angle in complete1.angles:
            angles_yes.add(angle)
        for angle in complete2.angles:
            angles_yes.add(angle)
        angles = angles_yes - angles_no
        for angle in angles:
            inst.angles.append(angle)

        dihedrals_yes = set()
        dihedrals_no = set()
        for dihedral in frag1.dihedrals:
            dihedrals_no.add(dihedral)
        for dihedral in frag2.dihedrals:
            dihedrals_no.add(dihedral)
        for dihedral in frag_prod.dihedrals:
            dihedrals_yes.add(dihedral)
        for dihedral in complete1.dihedrals:
            dihedrals_yes.add(dihedral)
        for dihedral in complete2.dihedrals:
            dihedrals_yes.add(dihedral)
        dihedrals = dihedrals_yes - dihedrals_no
        for dihedral in dihedrals:
            inst.dihedrals.append(dihedral)

        exclusions_yes = set()
        exclusions_no = set()
        for exclusion in frag1.exclusions:
            exclusions_no.add(exclusion)
        for exclusion in frag2.exclusions:
            exclusions_no.add(exclusion)
        for exclusion in frag_prod.exclusions:
            exclusions_yes.add(exclusion)
        for exclusion in complete1.exclusions:
            exclusions_yes.add(exclusion)
        for exclusion in complete2.exclusions:
            exclusions_yes.add(exclusion)
        exclusions = exclusions_yes - exclusions_no
        for exclusion in exclusions:
            inst.exclusions.append(exclusion)

        constraints_yes = set()
        constraints_no = set()
        for constraint in frag1.constraints:
            constraints_no.add(constraint)
        for constraint in frag2.constraints:
            constraints_no.add(constraint)
        for constraint in frag_prod.constraints:
            constraints_yes.add(constraint)
        for constraint in complete1.constraints:
            constraints_yes.add(constraint)
        for constraint in complete2.constraints:
            constraints_yes.add(constraint)
        constraints = constraints_yes - constraints_no
        for constraint in constraints:
            inst.constraints.append(constraint)

        impropers_yes = set()
        impropers_no = set()
        for improper in frag1.impropers:
            impropers_no.add(improper)
        for improper in frag2.impropers:
            impropers_no.add(improper)
        for improper in frag_prod.impropers:
            impropers_yes.add(improper)
        for improper in complete1.impropers:
            impropers_yes.add(improper)
        for improper in complete2.impropers:
            impropers_yes.add(improper)
        impropers = impropers_yes - impropers_no
        for improper in impropers:
            inst.impropers.append(improper)

    def modification(self, frag1, frag2, product):
        """Modification helper for the D/M algorithm
        """
        product_particles = frag1.particles + frag2.particles
        complete1 = self.destroy_fragment(frag1)
        complete2 = self.destroy_fragment(frag2)
        complete_particles = complete1.particles + complete2.particles
        assert len(complete_particles) >= len(product_particles)
        if len(complete_particles) > len(product_particles):
            frag_prod = self.instantiate_over_existing(
                product, product_particles, update=True
            )
            self.new_dynamic_complete_fragment(
                complete_particles, frag1, frag2, frag_prod, complete1,
                complete2
            )
        else:
            self.instantiate_over_existing(
                product, product_particles, first=True, update=True
            )

    def build_reaction_matrix(self):
        """Returns a hash table where reactions can be looked up for 2
        fragment names
        """
        return self.reactive_pairs

    def get_initiator_list(self):
        initiators: list[Fragment] = []
        for frag in self.frag_list:
            if frag is None:
                continue
            if frag.name in self.reactive_types:
                initiators.append(frag)

        return initiators

    def dump(self):
        print("==== TopStar / Fragment Types ====")
        for k, molfrag in self.type_lookup.items():
            sys.stdout.write(f"{k} ")
        sys.stdout.write("\n")
        print("==== TopStar / ReactionTemplates ====")
        for rx in self.reaction_list:
            print(rx)
        print("==== TopStar / Fragments ====")
        for frag in self.frag_list:
            print(frag)
