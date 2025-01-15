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
from .sysstar import SysStar
from .sysstar import Force, Interaction
from fnmatch import fnmatch
from .utils import pdist, backup_try

@dataclass
class FragFragment:
    name: str  # when instantiated it has this name
    mol: str
    name_id: str  # internally it has this name
    # list[(in_itp_id, type, name, is_edge)]
    atoms: list[(int, str, str, bool)]


@dataclass
class MolFragment:
    molecule_name: str
    # atoms: type, resnum, resname, atomname, chargegr, charge, mass
    atoms: list[(str, int, str, str, int, float, float)]
    # force i j params
    bonds: list[(Force, int, int, list)]
    # force i j k params
    angles: list[(Force, int, int, int, list)]
    # force i j k l params
    dihedrals: list[(Force, int, int, int, int, list)]
    # exclusions: list of id's
    exclusions: list[list[int]]
    # constraints: i j length
    constraints: list[(int, int, float)]
    # interactions: generic members and params
    interactions: list[(Force, list, list)]

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
                # I expect it to
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
    type: str  # one of: simple, poly, inter, intra, mono


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
    interactions: list[Interaction]
    exclusions: list[int]
    constraints: list[int]
    edge: list[bool]
    frag_id: int
    original_parent = None

    def __init__(self, name, id, original_parent=None):
        self.name = name
        self.particles = []
        self.interactions = []
        self.exclusions = []
        self.constraints = []
        self.edge = []
        self.frag_id = id
        self.original_parent = original_parent


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
        mol_fragment = MolFragment(name, [], [], [], [], [], [], [])
        self.type_lookup[name] = mol_fragment
        return mol_fragment

    def new_frag_fragment(self, name: str, parent: str):
        name_id = name
        i = 0
        while self.type_lookup.get(name_id):
            name_id = name + str(i)
            i += 1
        frag_fragment = FragFragment(name, parent, name_id, [])
        self.type_lookup[name_id] = frag_fragment
        if self.subfrag_map.get(parent):
            self.subfrag_map[parent].append(name_id)
        else:
            self.subfrag_map[parent] = [name_id]
        return frag_fragment

    def new_reaction(self, reaction: ReactionTemplate):
        self.reaction_list.append(reaction)
        self.reactive_pairs[(reaction.r1, reaction.r2)] = reaction
        self.reactive_pairs[(reaction.r2, reaction.r1)] = reaction
        self.reactive_types.add(reaction.r1)
        self.reactive_types.add(reaction.r2)

    def instantiate_subfrag(self, name_id, inst, original_parent=None):
        """Instantiates a subfragment subfrag, with forces in S* as seen in
        inst.

        args: subfrag: name_id, inst: Fragment"""
        frag: FragFragment = self.type_lookup.get(name_id)
        name = frag.name
        if frag is None:
            raise ValueError(f"Can't find frag {name_id}")
        elif not isinstance(frag, FragFragment):
            raise ValueError(f"Attempt to recursively instantiate {name_id}, "
                             "but it's not a frag fragment type."
                             f"It is: {name_id}")
        subinst_id = len(self.frag_list)
        subinst = Fragment(name, subinst_id, original_parent=original_parent)
        self.frag_list.append(subinst)
        normal_indices = set()
        all_indices = set()
        # particles
        for i, atom in enumerate(frag.atoms):
            # remapping of parent_id's to frag_id's
            parent_id, type, name, is_edge = atom
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
            if not is_edge:
                normal_indices.add(part_index)
            # build up subinst
            subinst.particles.append(part_index)
            subinst.edge.append(is_edge)
            # frag_id, in_frag_id, is_edge
            self.defrag_list[part_index].append((subinst_id, i, is_edge))

        # interactions get added if all participants are normal atoms
        for interaction in inst.interactions:
            members = interaction.get_members()
            include_interaction = True
            for member in members:
                if member not in normal_indices:
                    include_interaction = False
                    break
            if include_interaction:
                subinst.interactions.append(interaction)

        # exclusions, constaints get added if both of the atoms is a normal atom
        for excl_id in inst.exclusions:
            i, j = self.system.get_exclusion_members(excl_id)
            if i in normal_indices and j in normal_indices:
                subinst.exclusions.append(excl_id)
        for cid in inst.constraints:
            i, j = self.system.get_constraint_members(cid)
            if i in normal_indices and j in normal_indices:
                subinst.constraints.append(cid)

        # subfrags can contain further subfrags
        self.instantiate_subfrags(name, subinst)

    def instantiate_subfrags(self, frag_name, inst, original_parent=None):
        """Instantiates all subfrags of frag_name.
        inst is a Fragment instance that contains the reference to all forces
        in S*."""
        subfrags = self.subfrag_map.get(frag_name)
        if subfrags is not None:
            for subfrag in subfrags:
                self.instantiate_subfrag(subfrag, inst, original_parent)

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
            self, frag_name, particles, first=False, edge_list=False):
        """Takes a name of a mol fragment, creates new particles for it in
        the system and the corresponding interactions as well.
        Recursively instantiates all subfragments too.

        first: if set to true, it will insert into the defrag list, if false
        it will append

        edge_list: if false, won't change particle types, names, etc.
        if a list of booleans, it will change them where the booleans are False
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
            if edge_list:
                type, resnum, resname, atomname, chargegr, charge, mass = \
                    frag.atoms[in_frag_id]
                oldname, oldtype, oldcharge, oldmass = \
                    self.system.get_particle_details(part_id)
                # names are pattern matched rather than updated
                # so atom names actually stay the same as in monomers
                if not fnmatch(oldname, atomname):
                    raise Exception("Unmatching name during instantiate: "
                                    f"was {oldname}, pattern is {atomname}")
                # edge atoms are not updated
                # non edge atoms get updated, no match check for type, q, m
                if not edge_list[in_frag_id]:
                    self.system.update_particle(
                        part_id, oldname, type, charge, mass
                    )
        # bonds
        for (force, i, j, params) in frag.bonds:
            b = force.add(particles[i], particles[j], *params)
            inst.interactions.append(b)
        # angles
        for (force, i, j, k, params) in frag.angles:
            a = force.add(particles[i], particles[j], particles[k], *params)
            inst.interactions.append(a)
        # dihedrals
        for (force, i, j, k, l, params) in frag.dihedrals:
            d = force.add(
                particles[i], particles[j], particles[k], particles[l],
                *params
            )
            inst.interactions.append(d)
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
        # generic interactions
        for (force, members, params) in frag.interactions:
            members = [particles[x] for x in members]
            inst.interactions.append(force.add(members, params))

        self.instantiate_subfrags(frag_name, inst, original_parent=inst)
        return inst

    def get_molecule_fragment(self, frag: Fragment):
        # 0. save the molecule fragment
        molfrag_id, _, _ = self.defrag_list[frag.particles[0]][0]
        return self.frag_list[molfrag_id]

    def destroy_fragment(self, frag: Fragment):
        """args:
        frag: Fragment

        returns: the root molecule fragment from underneath
        """

        # 1. remove overlapping fragments
        # get the list of fragments to remove
        for part_id, part in enumerate(frag.particles):
            is_edge = frag.edge[part_id]
            defrag = self.defrag_list[part]
            i = 0
            while i < len(defrag):
                other_frag_id, _, is_other_edge = defrag[i]
                if not is_edge or not is_other_edge or \
                        other_frag_id == frag.frag_id:
                    # remove frag + remove other frags if the particle in frag
                    # is not edge OR if the other frag contains it as non edge
                    # 
                    # allowed to stay if it's edge-edge overlap
                    del defrag[i]
                    self.frag_list[other_frag_id] = None
                else:
                    i += 1

        # 2. remove forces in fragment
        for interaction in frag.interactions:
            interaction.remove()
        for excl in frag.exclusions:
            self.system.remove_exclusion(excl)
        for con in frag.constraints:
            # constraints can't be removed, so this will throw an exception
            self.system.remove_constraint(con)

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

        inter_yes = set()
        inter_no = set()
        for inter in frag1.interactions:
            inter_no.add(inter)
        for inter in frag2.interactions:
            inter_no.add(inter)
        for inter in frag_prod.interactions:
            inter_yes.add(inter)
        for inter in complete1.interactions:
            inter_yes.add(inter)
        for inter in complete2.interactions:
            inter_yes.add(inter)
        inters = inter_yes - inter_no
        for inter in inters:
            inst.interactions.append(inter)

        # TODO modularize exclusions and constraints into the same system later
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

    def detection(self,
                  frag1: Fragment,
                  frag2: Fragment,
                  rx: ReactionTemplate,
                  pos,
                  box
                  ):
        """Returns True if frag1 and frag2 fulfill constraints specified in rx

        Returns False if they should not react
        """
        # simple rules check
        if len(set(frag1.particles) & set(frag2.particles)) != 0:
            # Overlapping fragments can never react
            return False
        mol1 = self.get_molecule_fragment(frag1)
        mol2 = self.get_molecule_fragment(frag2)
        if rx.type == "simple":
            # simple bimolecular reaction with no additional restrictions
            pass
        elif rx.type == "poly":
            # subfragments instantiated together (same monomer)
            # cannot react with eachother
            if frag1.original_parent is not None and \
                    frag1.original_parent == frag2.original_parent:
                return False
        elif rx.type == "inter":
            # subfragments in the same molecule cannot react with
            # eachother
            if mol1 == mol2:
                return False
        elif rx.type == "intra":
            # subfragments only in the same molecule can react with
            # eachother
            if mol1 != mol2:
                return False
        elif rx.type == "mono":
            # there is only one fragment as reactant
            # frag1 == frag2 must be true i guess?
            # TODO this will take more work...
            if frag1 != frag2:
                return False
        else:
            raise ValueError(f"Unknown rx type {rx.type}")

        init1 = frag1.particles[0]
        init2 = frag2.particles[0]
        # position dependent checks
        dist = pdist(pos[init1], pos[init2], box)
        if dist >= rx.distance_max:
            return False
        return True

    def modification(self, frag1, frag2, product):
        """Modification helper for the D/M algorithm
        """
        # TODO check for overlapping frag1 and frag2 and reject such reactions
        complete1 = self.get_molecule_fragment(frag1)
        complete2 = self.get_molecule_fragment(frag2)
        if len(set(frag1.particles) & set(frag2.particles)) != 0:
            raise Exception("Overlapping fragments reacting")
        product_particles = frag1.particles + frag2.particles
        product_edge = frag1.edge + frag2.edge
        self.destroy_fragment(frag1)
        self.destroy_fragment(frag2)

        if complete1 != complete2:
            # bimolecular
            complete_particles = complete1.particles + complete2.particles
            assert len(complete_particles) >= len(product_particles)
            if len(complete_particles) > len(product_particles):
                frag_prod = self.instantiate_over_existing(
                    product, product_particles, edge_list=product_edge
                )
                self.new_dynamic_complete_fragment(
                    complete_particles, frag1, frag2, frag_prod, complete1,
                    complete2
                )
            else:
                self.instantiate_over_existing(
                    product, product_particles, first=True,
                    edge_list=product_edge
                )
        else:
            # intramolecular TODO
            frag_prod = self.instantiate_over_existing(
                product, product_particles, edge_list=product_edge
            )
            self.new_dynamic_complete_fragment(
                complete1.particles, frag1, frag2, frag_prod, complete1,
                complete1
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
        backup_try("top.dump")
        with open("top.dump", "w") as file:
            print("==== TopStar / Fragment Types ====", file=file)
            for k, molfrag in self.type_lookup.items():
                file.write(f"{k} ")
            file.write("\n")
            print("==== TopStar / ReactionTemplates ====", file=file)
            for rx in self.reaction_list:
                print(rx, file=file)
            print("==== TopStar / Fragments ====", file=file)
            for frag in self.frag_list:
                print(frag, file=file)
