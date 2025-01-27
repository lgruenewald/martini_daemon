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
from .forces.force import Force, Interaction
from fnmatch import fnmatch
from .utils import pdist, backup_try

@dataclass
class FragFragment:
    name: str  # when instantiated it has this name (e.g. in [ rx ])
    mol: str
    name_id: str  # internally it has this name for subfrag mapping
    # list[(in_itp_id, type, name, is_edge)]
    atoms: list[(int, str, str, bool)]


@dataclass
class MolFragment:
    molecule_name: str
    # atoms: type, resnum, resname, atomname, chargegr, charge, mass
    atoms: list[(str, int, str, str, int, float, float)]
    # TODO generalize this into just interactions
    # force i j params
    bonds: list[(Force, int, int, list)]
    # force i j k params
    angles: list[(Force, int, int, int, list)]
    # force i j k l params
    dihedrals: list[(Force, int, int, int, int, list)]
    # exclusions: list of id's
    exclusions: list[list[int]]
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
    break_pairs: list[(int, int)]
    update_groups: list[list[int]]


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
    frag_id: int

    def __init__(self, name, id):
        self.name = name
        self.particles = []
        self.interactions = []
        self.frag_id = id


class TopStar():
    # === Live molecule/fragment information ===
    # T* fragment and defrag list
    frag_list: dict[int, Fragment]
    next_frag_id: int
    # for every part_id have a list of fragments it is in
    # type: list[list[frag_id]]
    defrag_list: list[list[int]]

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
        self.frag_list = {}
        self.next_frag_id = 0
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
        mol_fragment = MolFragment(name, [], [], [], [], [], [])
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

    def add_frag_to_list(self, name: str):
        if name in self.reactive_types:
            inst = Fragment(name, self.next_frag_id)
            self.frag_list[self.next_frag_id] = inst
            self.next_frag_id += 1
            return inst
        else:
            return Fragment(name, -1)

    def instantiate_subfrag(self, name_id, inst):
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
        subinst = self.add_frag_to_list(name)
        all_indices = set()
        # particles
        for i, atom in enumerate(frag.atoms):
            # remapping of parent_id's to frag_id's
            parent_id, type, name = atom
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
            # build up these sets that are used for constructing the forces
            # prevent double inclusion of the same atom
            if part_index in all_indices:
                raise ValueError("Double inclusion of atom in subfrag")
            all_indices.add(part_index)
            # build up subinst
            subinst.particles.append(part_index)
            # frag_id, in_frag_id
            if subinst.frag_id != -1:
                self.defrag_list[part_index].append(subinst.frag_id)

        # interactions get added if all participants are normal atoms
        for interaction in inst.interactions:
            members = interaction.get_members()
            include_interaction = False
            for member in members:
                if member in all_indices:
                    include_interaction = True
                    break
            if include_interaction:
                subinst.interactions.append(interaction)

        # subfrags can contain further subfrags
        self.instantiate_subfrags(name, subinst)

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

    def instantiate_over_existing(self, frag_name, particles, reaction=False):
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
        inst = self.add_frag_to_list(frag_name)
        # particles
        for in_frag_id, part_id in enumerate(particles):
            inst.particles.append(part_id)
            # defrag
            if inst.frag_id != -1:
                self.defrag_list[part_id].append(inst.frag_id)
            if reaction:
                # update_types is only True if this is called during a reaction
                type, resnum, resname, atomname, chargegr, charge, mass = \
                    frag.atoms[in_frag_id]
                oldname, oldtype, oldcharge, oldmass = \
                    self.system.get_particle_details(part_id)
                # names are pattern matched rather than updated
                # so atom names actually stay the same as in monomers
                if not fnmatch(oldname, atomname):
                    raise Exception("Unmatching name during instantiate: "
                                    f"was {oldname}, pattern is {atomname}")
                # The * is only here as an option not to update types.
                if oldtype == type or type == "*":
                    # same type or type to remain same with *
                    # update if charge or mass change
                    if (charge is not None and oldcharge != charge) or \
                            (mass is not None and oldmass != mass):
                        self.system.update_particle(
                            part_id, oldname, oldtype, charge, mass
                        )
                else:
                    # new type, so gotta update anyway
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
                    e = self.system.exclusions.add(particles[i], particles[j])
                    inst.interactions.append(e)
        # generic interactions
        for (force, members, params) in frag.interactions:
            members = [particles[x] for x in members]
            inst.interactions.append(force.add(members, params))

        self.instantiate_subfrags(frag_name, inst)
        return inst

    def remove_fragment(self, frag: Fragment):
        """Removes a fragment from frag_lits and defrag_list
        """

        for part in frag.particles:
            defrag = self.defrag_list[part]
            i = 0
            while i < len(defrag):
                if defrag[i] == frag.frag_id:
                    del defrag[i]
                else:
                    i += 1

        del self.frag_list[frag.frag_id]

    def destroy_fragment(self, frag: Fragment):
        """args:
        frag: Fragment

        removes overlapping fragments, then removes fragment
        """

        for part in frag.particles:
            # "cache", since remove_fragment will change this list
            defrag = self.defrag_list[part].copy()
            for other_frag_id in defrag:
                # remove self at the end, not here
                if frag.frag_id == other_frag_id:
                    continue
                # if the other frag hasn't been removed yet, remove it
                other_frag = self.frag_list.get(other_frag_id)
                if other_frag is not None:
                    self.remove_fragment(other_frag)

        self.remove_fragment(frag)

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

        init1 = frag1.particles[0]
        init2 = frag2.particles[0]
        # position dependent checks
        dist = pdist(pos[init1], pos[init2], box)
        if dist > rx.distance_max:
            return False
        return True

    def modification(self, frag1, frag2, product):
        """Modification helper for the D/M algorithm
        """
        product_particles = frag1.particles + frag2.particles
        self.destroy_fragment(frag1)
        self.destroy_fragment(frag2)

        self.instantiate_over_existing(
            product, product_particles, reaction=True
        )

    def build_reaction_matrix(self):
        """Returns a hash table where reactions can be looked up for 2
        fragment names
        """
        return self.reactive_pairs

    def get_initiator_list(self):
        initiators: list[Fragment] = []
        for id, frag in self.frag_list.items():
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
            for id, frag in self.frag_list.items():
                print(id, frag, file=file)
