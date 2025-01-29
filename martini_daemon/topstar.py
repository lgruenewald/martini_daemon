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
from .utils import pdist, backup_try, pcos_angle, pdihedral
import random
import numpy as np


random.seed()


@dataclass
class FragFragment:
    name: str  # when instantiated it has this name (e.g. in [ rx ])
    mol: str
    name_id: str  # internally it has this name for subfrag mapping
    # list[(in_itp_id, type, name)]
    atoms: list[(int, str, str)]


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


class ReactionTemplate:
    name: str
    r1: str
    r2: str
    p1: str
    distance_max: list[(int, int, float)]
    distance_min: list[(int, int, float)]
    probability: float
    global_counter: int
    global_limit: int
    break_groups: list[list[int]]
    update_groups: list[list[int]]

    def __init__(self, name):
        self.name = name
        self.r1 = None
        self.r2 = None
        self.p1 = None
        self.distance_max = []
        self.distance_min = []
        self.angle_limits = []
        self.dihedral_limits = []
        self.probability = 1.0
        self.global_counter = 0
        self.global_limit = None
        self.break_groups = []
        self.update_groups = []

    def is_complete(self):
        # when a reaction is finished parsing, if this returns False
        # it is considered an error

        # currently: must have a reactant and product and at least one distance
        # constraint, but this is a bit arbitrary / should reflect common
        # potential mistakes people would make when writing *.rx files
        return (
            self.r1 is not None and
            self.p1 is not None and
            (len(self.distance_max) > 0 or len(self.distance_min) > 0)
        )


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

    log_path: str

    def __init__(self, system, log_path=None):
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
        self.log_path = log_path
        if log_path is not None:
            backup_try(log_path)

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
# TODO: find out why this inverted thing worked in the past
#        self.reactive_pairs[(reaction.r2, reaction.r1)] = reaction
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

        self.log(f"[info] instantiate_subfrag {frag.name} over particles {frag.atoms}")
        # interactions get added if all participants are normal atoms
        # TODO functionalize these checks
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

    def instantiate_over_existing(
            self, frag_name, particles, reaction=False,
            interactions=[]):
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
        if reaction:
            self.log("instantiate_over_existing during reaction "
                     f"product_name {frag_name} "
                     f"transferring interactions #: {len(interactions)}")

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
        # interactions inherited / not added here
        for inter in interactions:
            inst.interactions.append(inter)

        self.instantiate_subfrags(frag_name, inst)
        return inst

    def remove_fragment(self, frag: Fragment):
        """Removes a fragment from frag_lits and defrag_list
        """

        # TODO functionalize these checks
        for part in frag.particles:
            defrag = self.defrag_list[part]
            i = 0
            while i < len(defrag):
                if defrag[i] == frag.frag_id:
                    del defrag[i]
                else:
                    i += 1

        del self.frag_list[frag.frag_id]

    def destroy_fragment(self, frags: list[Fragment], modified_atoms: list[int]):
        """args:
        frag: Fragment

        removes overlapping fragments with modified_atoms,
        then removes fragments in frags
        """

        frag_ids = {frag.frag_id for frag in frags}

        # TODO functionalize these checks
        for part in modified_atoms:
            # "cache", since remove_fragment will change this list
            defrag = self.defrag_list[part].copy()
            for other_frag_id in defrag:
                # remove self at the end, not here
                if other_frag_id in frag_ids:
                    continue
                # if the other frag hasn't been removed yet, remove it
                other_frag = self.frag_list.get(other_frag_id)
                if other_frag is not None:
                    self.remove_fragment(other_frag)

        for frag in frags:
            self.remove_fragment(frag)

    def pre_detection(self):
        self.log("[info] pre_detection hook")
        for rx in self.reaction_list:
            rx.global_counter = 0

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
        # the order of checks should be from fastest to slowest to maximize
        # performance
        # simple rules check
        # fragments that were removed cannot react any more

        # overlapping fragments can never react:
        if frag2 is not None:
            for i in frag1.particles:
                for j in frag2.particles:
                    if i == j:
                        return False

        # limiter checks, first for performance
        if rx.global_limit is not None and rx.global_counter >= rx.global_limit:
            return False
        if random.random() > rx.probability:
            return False

        # position dependent checks
        particles = frag1.particles.copy()
        if frag2 is not None:
            particles += frag2.particles

        for (i, j, rmax) in rx.distance_max:
            init1 = particles[i]
            init2 = particles[j]
            dist = pdist(pos[init1], pos[init2], box)
            if dist > rmax:
                return False

        for (i, j, rmin) in rx.distance_min:
            init1 = particles[i]
            init2 = particles[j]
            dist = pdist(pos[init1], pos[init2], box)
            if dist < rmin:
                return False

        for (i, j, k, cos_min, cos_max) in rx.angle_limits:
            # the particle positions of particle i, j, k
            p1 = pos[particles[i]]
            p2 = pos[particles[j]]
            p3 = pos[particles[k]]
            cos = pcos_angle(p1, p2, p3, box)
            if cos > cos_min and cos < cos_max:
                self.log(f"Rejected reaction because of cos_angle {cos} "
                         f"inside of {cos_min} to {cos_max} range")
                return False

        for (i, j, k, l, min, max) in rx.dihedral_limits:
            # particle positions
            p1 = pos[particles[i]]
            p2 = pos[particles[j]]
            p3 = pos[particles[k]]
            p4 = pos[particles[l]]
            theta = pdihedral(p1, p2, p3, p4, box)
            if theta > min and theta < max:
                self.log(f"Rejected reaction because of dihedral {theta} "
                         f"inside of {min} to {max} range")
                return False

        rx.global_counter += 1
        return True

    def pre_modification(self):
        # hook that gets called after detection, before modification
        # only called if there is any modification going on
        self.log("[info] pre modification hook")

    def modification(self, frag1, frag2, rx: ReactionTemplate):
        """Modification helper for the D/M algorithm
        """
        self.log("[info] Modification algo!")
        self.log(f"reaction {rx.name}")
        product = rx.p1
        product_particles = frag1.particles.copy()
        all_interactions = frag1.interactions.copy()
        self.log(f"frag1 particles {frag1.particles}")
        if frag2 is not None:
            self.log(f"frag2 particles {frag2.particles}")
            product_particles += frag2.particles
            all_interactions += frag2.interactions
        self.log(f"product_particles {product_particles}")
        modified_atoms = set()
        i = 0
        while i < len(all_interactions):
            # TODO functionalize these checks
            interaction = all_interactions[i]
            remove = False
            members = interaction.get_members()
            # process rx_break
            # if any in group not in members -> not remove candidate
            for group in rx.break_groups:
                all = True
                for atom in group:
                    if atom not in members:
                        all = False
                        break
                if all:
                    remove = True

            # process rx_update
            # if any in members not in group -> not remove candidate
            for group in rx.update_groups:
                all = True
                for atom in members:
                    if atom not in group:
                        all = False
                        break
                if all:
                    remove = True
            if remove:
                # mark atoms that were modified for fragment overlap deleting
                for atom in members:
                    modified_atoms.add(atom)
                interaction.remove()
                del all_interactions[i]
            else:
                i += 1

        modified_atoms = list(modified_atoms)
        self.log(f"modified_atoms {modified_atoms}")

        if frag2 is not None:
            self.destroy_fragment([frag1, frag2], modified_atoms)
        else:
            self.destroy_fragment([frag1], modified_atoms)

        self.instantiate_over_existing(
            product, product_particles, reaction=True,
            interactions=all_interactions
        )

    def post_modification(self):
        # hook that only gets called after modification
        self.log("[info] post modification hook")
        if self.log_path is not None:
            self.dump(self.log_path, True)

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

    def log(self, message):
        if self.log_path is not None:
            with open(self.log_path, "a") as file:
                print(message, file=file)

    def dump(self, path="top.dump", append=False):
        if not append:
            backup_try(path)
        with open(path, "a") as file:
            print("===== TopStar Dump =====", file=file)
            print("==== TopStar / Fragment Types ====", file=file)
            for k, molfrag in self.type_lookup.items():
                file.write(f"{k} ")
            file.write("\n")
            print("==== TopStar / ReactionTemplates ====", file=file)
            for rx in self.reaction_list:
                print(f"rx {rx.name} r1 {rx.r1} r2 {rx.r2} p1 {rx.p1}", file=file)
            print("==== TopStar / Fragments ====", file=file)
            for id, frag in self.frag_list.items():
                print(f"{id}: <frag {frag.name} ps {frag.particles}>", file=file)
