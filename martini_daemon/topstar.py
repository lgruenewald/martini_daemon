#!/usr/bin/env python3
from dataclasses import dataclass
import random
from freud.box import Box
from freud.locality import NeighborList, AABBQuery
import numpy as np
from typing import Optional
from .sysstar import SysStar
from .forces.force import Force, Interaction
from .utils import pdist, pcos_angle, pdihedral
from .reporters.reporter import Reporter
from .graph import GraphFragment, GraphMatch, GraphAtomType
from .fragment import Fragment, Fragments

random.seed()


@dataclass
class NumberedFragment:
    name: str  # when instantiated it has this name (e.g. in [ rx ])
    mol: str
    name_id: str  # internally it has this name for subfrag mapping
    # list[in_itp_id]
    atoms: list[int]


# TODO move this to its own file and move all the per fragment type helpers
# (e.g.) instantiate to their files, break up T* monolithicness a bit...
class MolFragment:
    molecule_name: str
    # atoms: type, resnum, resname, atomname, chargegr, charge, mass
    atoms: list[tuple[str, int, str, str, int, float, float]]
    # tuple[int, str] instead of int when reaction TODO
    exclusions: set[tuple[int, int]]
    # interactions: generic members and params
    # tuple[int, str] instead of int when reaction
    interactions: list[tuple[Force, list[int], list[float]]]
    index_type = "index"

    def __init__(self, name):
        self.molecule_name = name
        self.atoms = []
        self.exclusions = set()
        self.interactions = []
        self.renames = []
        self.retypes = []

    def add_exclusion(self, i, j) -> None:
        """Adds an exclusion to the list of exclusions
        Does not add duplicate exclusions."""

        self.exclusions.add((i, j))
        self.exclusions.add((j, i))


class ReactionTemplate:
    name: str
    reactants: list[str]
    products: list[tuple[str, Optional[list[tuple[int, str]]]]]
    distance_max: list[tuple[tuple[int, str], tuple[int, str], float]]
    distance_min: list[tuple[tuple[int, str], tuple[int, str], float]]
    angle_limits: list[tuple[tuple[int, str], tuple[int, str], tuple[int, str],
                       float, float]]
    dihedral_limits: list[tuple[tuple[int, str], tuple[int, str],
                          tuple[int, str], tuple[int, str], float, float]]
    probability: float
    global_counter: int
    global_limit: int
    break_groups: list[list[tuple[int, str]]]
    update_groups: list[list[tuple[int, str]]]
    skip: int
    product: MolFragment
    renames: list[tuple[int, str, str]]
    retypes: list[tuple[int, str, str]]
    recharges: list[tuple[int, str, float]]
    remasses: list[tuple[int, str, float]]

    def __init__(self, name):
        self.name = name
        self.reactants = []
        self.products = []
        self.distance_max = []
        self.distance_min = []
        self.angle_limits = []
        self.dihedral_limits = []
        self.probability = 1.0
        self.global_counter = 0
        self.global_limit = None
        self.break_groups = []
        self.update_groups = []
        self.skip = None
        self.product = MolFragment(name)
        self.product.index_type = "pair"
        self.renames = []
        self.retypes = []
        self.recharges = []
        self.remasses = []

    def is_complete(self) -> bool:
        # when a reaction is finished parsing, if this returns False
        # it is considered an error

        return (
            len(self.reactants) > 0
        )


class TopStar():
    # ======= (1/3) Building things =======
    # T* fragment and defrag list
    frag_list: dict[int, Fragment]
    next_frag_id: int
    # for every part_id have a list of fragments it is in
    # and a list of interactions it is in
    # type is  list[list[frag_id]]
    defrag_list: list[list[int]]
    interaction_list: list[list[Interaction]]

    # type name -> list of subfrag names
    subfrag_map: dict[str, list[str]]

    graph_fragment_list: list[GraphFragment]

    # type name -> type, used for instantiation
    # TODO only use it for MolFragment
    type_lookup: dict[str, NumberedFragment | MolFragment]

    reaction_list: list[ReactionTemplate]
    reactive_types: set[str]

    system: SysStar

    reporters: list[Reporter]

    def __init__(self, system, logger):
        self.frag_list = {}
        self.next_frag_id = 0
        self.defrag_list = []
        self.interaction_list = []
        self.frag_fragments = []
        self.mol_fragments = []
        self.reaction_list = []
        self.type_lookup = {}
        self.subfrag_map = {}
        self.reactive_types = set()
        self.reporters = []
        self.graph_fragment_list = []
        self.system = system
        self.logger = logger

    def add_reporter(self, reporter) -> None:
        self.reporters.append(reporter)

    def new_mol_fragment(self, name: str) -> MolFragment:
        if self.type_lookup.get(name):
            raise ValueError(f"Second definition of fragment type {name}")
        mol_fragment = MolFragment(name)
        self.type_lookup[name] = mol_fragment
        return mol_fragment

    def new_numbered_fragment(self, name: str, parent: str,
                              particles: list[int]) -> None:
        name_id = name
        i = 0
        while self.type_lookup.get(name_id):
            name_id = name + str(i)
            i += 1
        frag = NumberedFragment(name, parent, name_id, particles)
        self.type_lookup[name_id] = frag
        if self.subfrag_map.get(parent):
            self.subfrag_map[parent].append(name_id)
        else:
            self.subfrag_map[parent] = [name_id]

    def new_graph_fragment(self, frag: GraphFragment) -> None:
        self.graph_fragment_list.append(frag)

    def new_reaction(self, reaction: ReactionTemplate) -> MolFragment:
        self.reaction_list.append(reaction)
        for r in reaction.reactants:
            self.reactive_types.add(r)
        if self.type_lookup.get(reaction.name):
            raise ValueError(f"Second definition of fragment type / reaction {name}")
        self.type_lookup[reaction.name] = reaction.product
        return reaction.product

    def add_frag_to_list(self, name: str) -> Fragment:
        if name in self.reactive_types:
            inst = Fragment(name, self.next_frag_id)
            self.frag_list[self.next_frag_id] = inst
            self.next_frag_id += 1
            return inst
        else:
            return Fragment(name, -1)

    # Graph helpers

    def try_match_graphs(self, particles: set[int], molname=None) -> None:
        """
            Graph version of instantiate subfrag.

            Tries to match all known graphs at all particles in a list.
            Should be called after interactions (self.interaction_list) have
            been updated. Particles should be all affected particles, and
            all their neighbors (during reactions) or all particles
            (at the start).

            Will add graphs to self.frag_list, if the same match/graph doesn't
            already exist.

            If molname is specified, it means that we can be assured that
            the same matches are going to happen when we call it with the
            same molname again, so we can cache the results and speed up
            this function call later.
        """

        matches: list[GraphMatch] = []
        for graph in self.graph_fragment_list:
            # for each possible graph to match
            # don't match if it's only the specific molecule
            if len(graph.molecules) > 0:
                if molname not in graph.molecules:
                    continue
            matches += graph.match_particles(
                particles, self.system,
                self.interaction_list
            )

        for m in matches:
            inst = self.add_frag_to_list(m.graph.name)
            inst.graph = m.graph
            for key, _, _, type in m.graph.atoms:
                val = m.atoms.get(key)
                # key - name in the graph
                # val - particle id
                if val is None:
                    if type == GraphAtomType.OPT:
                        inst.opt[key] = None
                    continue
                if type == GraphAtomType.NORMAL:
                    inst.particles[key] = val
                elif type == GraphAtomType.OPT:
                    inst.opt[key] = val
                else:
                    raise NotImplementedError
                if inst.frag_id != -1:
                    self.defrag_list[val].append(inst.frag_id)

    def instantiate_subfrag(self, name_id: str, inst: Fragment) -> None:
        """Instantiates a subfragment subfrag, with forces in S* as seen in
        inst.

        args: subfrag: name_id, inst: Fragment"""
        frag: NumberedFragment = self.type_lookup.get(name_id)
        name = frag.name
        if frag is None:
            raise ValueError(f"Can't find frag {name_id}")
        elif not isinstance(frag, NumberedFragment):
            raise ValueError(f"Attempt to recursively instantiate {name_id}, "
                             "but it's not a frag fragment type."
                             f"It is: {name_id}")
        subinst = self.add_frag_to_list(name)
        all_indices = set()
        # particles
        for i, parent_id in enumerate(frag.atoms):
            # remapping of parent_id's to frag_id's
            # parent_id must be i range
            if parent_id < 0 or parent_id >= len(inst.particles):
                raise ValueError("Parent particle ID out of range")
            part_index = inst.index_atom(parent_id)
            # build up these sets that are used for constructing the forces
            # prevent double inclusion of the same atom
            if part_index in all_indices:
                raise ValueError("Double inclusion of atom in subfrag")
            all_indices.add(part_index)
            # build up subinst
            subinst.particles[f"{i+1}"] = part_index
            # frag_id, in_frag_id
            if subinst.frag_id != -1:
                self.defrag_list[part_index].append(subinst.frag_id)

        # subfrags can contain further subfrags
        self.instantiate_subfrags(name, subinst)

    def instantiate_subfrags(self, frag_name, inst) -> None:
        """Instantiates all subfrags of frag_name.
        inst is a Fragment instance that contains the reference to all forces
        in S*."""
        subfrags = self.subfrag_map.get(frag_name)
        if subfrags is not None:
            for subfrag in subfrags:
                self.instantiate_subfrag(subfrag, inst)

    def instantiate(self, frag_name: str) -> Fragment:
        # frag_name must refer to a MolFragment type
        frag = self.type_lookup.get(frag_name)
        if frag is None:
            raise ValueError(f"Can't find mol {frag_name}")
        elif not isinstance(frag, MolFragment):
            raise ValueError(f"Attempt to instantiate {frag_name}, but it's"
                             " not a mol fragment type."
                             f" It is: {frag}")
        # add particles to S*
        parts = []
        prev_resnum = 0
        for in_frag_id, atom in enumerate(frag.atoms):
            type, resnum, resname, atomname, chargegr, charge, mass = atom
            if resnum != prev_resnum:
                self.system.new_residue()
                prev_resnum = resnum
            p = self.system.add_particle(atomname, resname, type, charge, mass)
            parts.append(p)
            self.defrag_list.append([])
            self.interaction_list.append([])
        # instantiate interactions
        res = self.instantiate_over_existing(frag,
                                             Fragments.from_parts(parts))
        # add graphs to system
        self.try_match_graphs(set(parts), frag_name)
        return res

    def instantiate_over_existing(self, molfrag: MolFragment,
                                  frags: Fragments,
                                  reaction: bool = False,
                                  template: Optional[ReactionTemplate] = None
                                  ) -> Fragment:
        """Takes a name of a mol fragment, adds interactions to those particles
        according to the mol fragment.
        Recursively instantiates all subfragments too.
        """

        inst = self.add_frag_to_list(molfrag.molecule_name)
        # particles
        for in_frag_id, part_id in enumerate(frags.particles()):
            inst.particles[f"{in_frag_id+1}"] = part_id
            # defrag
            if inst.frag_id != -1:
                self.defrag_list[part_id].append(inst.frag_id)
            if reaction:
                # old system:
                # products also change type, charge, mass but not name
                type, _, _, _, _, charge, mass = \
                    molfrag.atoms[in_frag_id]
                # Atom names are not updated.
                # The * is only here as an option not to update types.
                if type != "*":
                    self.system.retype(part_id, type)
                if charge is not None:
                    self.system.recharge(part_id, charge)
                if mass is not None:
                    self.system.remass(part_id, mass)
        # new system: re* overrides everything though
        if template is not None:
            for (pair, new_name) in template.renames:
                found, part_id = frags.get(pair)
                if not found:
                    raise ValueError(f"Index out of range {pair}")
                if part_id is None:
                    continue
                self.system.rename(part_id, new_name)
            for (pair, new_type) in template.retypes:
                found, part_id = frags.get(pair)
                if not found:
                    raise ValueError(f"Index out of range {pair}")
                if part_id is None:
                    continue
                self.system.retype(part_id, new_type)
            for (pair, new_charge) in template.recharges:
                found, part_id = frags.get(pair)
                if not found:
                    raise ValueError(f"Index out of range {pair}")
                if part_id is None:
                    continue
                self.system.recharge(part_id, new_charge)
            for (pair, new_mass) in template.remasses:
                found, part_id = frags.get(pair)
                if not found:
                    raise ValueError(f"Index out of range {pair}")
                if part_id is None:
                    continue
                self.system.remass(part_id, new_charge)
        # exclusions
        for (i, j) in molfrag.exclusions:
            foundi, pi = frags.get(i)
            foundj, pj = frags.get(j)
            if not foundi or not foundj:
                # wrong atom name
                raise ValueError("Exclusion index out of range")
            if pi is None or pj is None:
                # optional atom missing
                continue
            if pi < pj:
                e = self.system.exclusions.add(pi, pj)
                self.interaction_list[pi].append(e)
                self.interaction_list[pj].append(e)
        # generic interactions
        for (force, members, params) in molfrag.interactions:
            member_parts = []
            for x in members:
                found, part = frags.get(x)
                if not found:
                    raise ValueError("Wrong atom name")
                member_parts.append(part)
            if any(map(lambda x: x is None, member_parts)):
                # optional missing => skip
                continue
            f = force.add(member_parts, params)
            for member in member_parts:
                self.interaction_list[member].append(f)

        self.instantiate_subfrags(molfrag.molecule_name, inst)
        return inst

    # ======= (2/3) Detection things =======
    def pre_detection(self, i):
        for rx in self.reaction_list:
            rx.global_counter = 0
        for reporter in self.reporters:
            reporter.pre_detection(i)

    def detection(self,
                  reactants: list[Fragment],
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
        pset = set()
        for r in reactants:
            for p in r.particles.values():
                if p in pset:
                    return False
                pset.add(p)
            for p in filter(lambda x: x is not None, r.opt.values()):
                if p in pset:
                    return False
                pset.add(p)

        # limiter checks, first for performance
        if rx.global_limit is not None and rx.global_counter >= rx.global_limit:
            return False
        if random.random() > rx.probability:
            return False

        # position dependent checks
        for ((moli, idi), (molj, idj), rmax) in rx.distance_max:
            found1, init1 = reactants[moli].get_atom(idi)
            found2, init2 = reactants[molj].get_atom(idj)
            if not found1 or not found2:
                # wrong atom names
                raise ValueError("Bad atom name in reaction condition")
            if init1 is None or init2 is None:
                # missing optional atoms
                continue
            dist = pdist(pos[init1], pos[init2], box)
            if dist > rmax:
                return False

        for ((moli, idi), (molj, idj), rmin) in rx.distance_min:
            found1, init1 = reactants[moli].get_atom(idi)
            found2, init2 = reactants[molj].get_atom(idj)
            if not found1 or not found2:
                # wrong atom names
                raise ValueError("Bad atom name in reaction condition")
            if init1 is None or init2 is None:
                # missing optional atoms
                continue
            dist = pdist(pos[init1], pos[init2], box)
            if dist < rmin:
                return False

        for ((moli, idi), (molj, idj), (molk, idk), cos_min, cos_max) in \
                rx.angle_limits:
            found1, p1 = reactants[moli].get_atom(idi)
            found2, p2 = reactants[molj].get_atom(idj)
            found3, p3 = reactants[molk].get_atom(idk)
            if not found1 or not found2 or not found3:
                # wrong atom names
                raise ValueError("Bad atom name in reaction condition")
            if p1 is None or p2 is None or p3 is None:
                # missing optional atoms
                continue
            # the particle positions of particle i, j, k
            cos = pcos_angle(pos[p1], pos[p2], pos[p3], box)
            if cos < cos_min and cos > cos_max:
                # inverted comparison because cosine is a constantly decreasing
                # function, cos_min is the minimum angle => max cosine value
                # cos_max is the maximum angle => min cosine value
                return False

        for ((moli, idi), (molj, idj), (molk, idk), (moll, idl), min, max) in \
                rx.dihedral_limits:
            # particle positions
            found1, p1 = reactants[moli].get_atom(idi)
            found2, p2 = reactants[molj].get_atom(idj)
            found3, p3 = reactants[molk].get_atom(idk)
            found4, p4 = reactants[moll].get_atom(idl)
            if not found1 or not found2 or not found3 or not found4:
                # wrong atom names
                raise ValueError("Bad atom name in reaction condition")
            if p1 is None or p2 is None or p3 is None or p4 is None:
                # missing optional atoms
                continue
            theta = pdihedral(pos[p1], pos[p2], pos[p3], pos[p4], box)
            if theta > min and theta < max:
                return False

        rx.global_counter += 1
        return True

    def get_init_map(self):
        init_map: dict[str, list[int]] = {}
        for key, frag in self.frag_list.items():
            if init_map.get(frag.name) is None:
                init_map[frag.name] = []
            init_map[frag.name].append(key)

        return init_map

    def get_neighbor_list(self, box, pos) -> \
            tuple[NeighborList, dict[int, int], dict[int, int]]:
        if len(self.frag_list) == 0:
            return None, {}, {}
        pos_filtered = []
        frag_to_nlist = {}
        nlist_to_frag = {}
        for id, frag in self.frag_list.items():
            filtered_id = len(pos_filtered)
            frag_to_nlist[id] = filtered_id
            nlist_to_frag[filtered_id] = id
            part_id = frag.index_atom(0)
            pos_filtered.append(pos[part_id])
        pos_filtered = np.array(pos_filtered)
        box = Box(box[0], box[1], box[2])
        query = NeighborList()
        query = AABBQuery(box, pos_filtered)
        return query, frag_to_nlist, nlist_to_frag

    def get_reaction_tree(self) -> dict:
        tree = {}
        for rx in self.reaction_list:
            reactants = rx.reactants
            node = tree
            for j, reactant in enumerate(reactants):
                if rx.skip is not None and j >= rx.skip:
                    # unelegant way of adding enum-like info in the k/v mapping
                    # ideally this would be possible at the type level
                    reactant = "*" + reactant
                if node.get(reactant) is None:
                    node[reactant] = {}
                node = node[reactant]
            node["_reaction"] = rx
        return tree

    # TODO type this better
    def detection_over_types(
            self,
            reactions: list[tuple[list[Fragment], ReactionTemplate]],
            skip: set[int],
            node: dict,
            previous_types: list[str],
            previous_indices: list[int],
            init_map, pos, box, query: AABBQuery,
            frag_to_nlist: dict[int, int], nlist_to_frag: dict[int, int]
            ) -> tuple[bool, list, set]:
        """
            Runs the detection algorithm over a single node in the reaction
            "tree".

            reactions - dynamic list of reactions that passed the D algo
            skip - dynamic set of indices in self.frag_list that already reacted
            so they cannot any more
            node - (sub)tree of reaction types left to check
            previous_types - reactant types already checked
            previous_indicies - indicies in init_map[reactant_type] that were
            already checked
            init_map - dict of reactant_type -> list of indicies of that type
            in self.frag_list
            pos - numpy array of positions
            box - periodic box info
        """

        for next_reactant, next_node in node.items():
            if next_reactant == "_reaction":
                # convert init_map indices to init_list indices
                frag_ids = [init_map[name.lstrip("*")][id] for name, id in
                            zip(previous_types, previous_indices)]
                frags = [self.frag_list[i] for i in frag_ids]
                rx = next_node
                if self.detection(frags, rx, pos, box):
                    # only until rx.skip do we add things to skip and rxs
                    # the rest are only for self.detection but nothing after
                    if rx.skip is not None:
                        frag_ids = frag_ids[:rx.skip]
                        frags = frags[:rx.skip]
                    for frag_id in frag_ids:
                        skip.add(frag_id)
                    reactions.append((frags, rx))
                    return reactions, skip, True

            if init_map.get(next_reactant.lstrip("*")) is None:
                continue

            # find out where to start indexing from
            # this also works well with the * hack for rx.skip
            # this used to be the if j <= i: continue check
            # prevents both self reactions and double counting
            start_at = 0
            for i, prev_reactant in enumerate(previous_types):
                if next_reactant == prev_reactant:
                    start_at = max(start_at, previous_indices[i] + 1)
            # go over all options for the next_reactant type
            # TODO rewrite this function completely with a more elegant
            # neighbor list
            neighbors = None
            if len(previous_indices) > 0:
                last_type = previous_types[-1]
                last_index = init_map[last_type][previous_indices[-1]]
                last_frag = self.frag_list[last_index]
                last_part = last_frag.index_atom(0)
                last_pos = np.array([pos[last_part]])
                nlist = query.query(last_pos, {"r_max": 3.}).toNeighborList()
                neighbors = set()
                for _, j in nlist[:]:
                    neighbors.add(nlist_to_frag[j])
            for j in range(start_at, len(init_map[next_reactant.lstrip("*")])):
                # skip is a set of init_list indices that already reacted
                # only check it if the current reactant is before rx.skip => not *
                if next_reactant[0] != "*" and init_map[next_reactant][j] in skip:
                    continue
                if neighbors is not None and j not in neighbors:
                    continue
                previous_types.append(next_reactant)
                previous_indices.append(j)
                reactions, skip, reacted = self.detection_over_types(
                    reactions, skip, next_node,
                    previous_types, previous_indices,
                    init_map, pos, box, query, frag_to_nlist, nlist_to_frag
                )
                previous_types.pop(-1)
                previous_indices.pop(-1)
                if reacted and len(previous_indices) > 0:
                    # there was a reaction, so previous_indices[0] is now
                    # skipped, so we can jump up all the way to root
                    return reactions, skip, True

        return reactions, skip, False

    # ======= (3/3) Modification things =======
    def clean_defrag(self, frag: Fragment, part: int) -> None:
        self.defrag_list[part] = list(filter(
            lambda x: x != frag.frag_id,
            self.defrag_list[part]
        ))

    def remove_fragment(self, frag: Fragment) -> None:
        """Removes a fragment from frag_list and defrag_list
        """

        for part in frag.particles.values():
            self.clean_defrag(frag, part)

        for part in filter(lambda x: x is not None, frag.opt.values()):
            self.clean_defrag(frag, part)

        del self.frag_list[frag.frag_id]

    def remove_interaction(self, interaction: Interaction) -> None:
        """Removes an interaction from interaction_list and S*
        """
        for part in interaction.get_members():
            self.interaction_list[part] = list(filter(
                lambda x: x != interaction,
                self.interaction_list[part]
            ))
        interaction.remove()

    def process_break(self, frags: list[Fragment], rx: ReactionTemplate,
                      ) -> None:
        """
            Process [rx_break] in rx over frags.
        """

        for group in rx.break_groups:
            group_atoms = []
            all_found = True
            for moli, idi in group:
                found, part = frags[moli].get_atom(idi)
                if not found:
                    raise ValueError("Unknown part name")
                if part is None:
                    all_found = False
                    break
                group_atoms.append(part)
            if not all_found:
                # missing optional
                continue
            # breaking interactions
            # for every group in break_groups
            # check all interactions and break if all in group are in inter
            for inter in self.interaction_list[group_atoms[0]][:]:
                if all(map(
                           lambda group_member:
                           group_member in inter.get_members(),
                           group_atoms
                       )):
                    self.remove_interaction(inter)

    def process_update(self, frags: list[Fragment], rx: ReactionTemplate,
                       ) -> None:
        """
            Process [rx_update] in rx over frags.
        """
        # TODO: consider if we really need this or if there are better ways
        # to remove interactions

        for group in rx.update_groups:
            group_atoms = []
            all_found = True
            for moli, idi in group:
                found, part = frags[moli].get_atom(idi)
                if not found:
                    raise ValueError("Unknown part name")
                if part is None:
                    all_found = False
                    break
                group_atoms.append(part)
            if not all_found:
                # missing optional
                continue
            # breaking interactions
            # check that all in interaction are in group
            for inter in self.interaction_list[group_atoms[0]][:]:
                if all(map(
                           lambda inter_member: inter_member in group_atoms,
                           inter.get_members()
                       )):
                    self.remove_interaction(inter)

    def populate_neighbors(self, particles: set[int]) -> set[int]:
        res = set()
        for part in particles:
            res.add(part)
            for inter in self.interaction_list[part]:
                for member in inter.get_members():
                    res.add(member)
        return res

    def remove_overlapping_graphs(self, particles: set[int]) -> None:
        for part in particles:
            for frag_id in self.defrag_list[part][:]:
                frag = self.frag_list.get(frag_id)
                if frag is not None and frag.graph is not None:
                    self.remove_fragment(frag)

    def pre_modification(self, rx_list: list[(list, ReactionTemplate)], i
                         ) -> None:
        # hook that gets called after detection, before modification
        # only called if there is any modification going on
        for reporter in self.reporters:
            reporter.pre_modification(rx_list, i)

    def modification(self, frags: list[Fragment], rx: ReactionTemplate
                     ) -> None:
        """Modification helper for the D/M algorithm
        """
        # just normal particles to instantiate products over
        product_particles = []
        # which particles to recalculate graphs over
        graph_recalc = set()
        for f in frags:
            product_particles += list(f.particles.values())
            graph_recalc |= set(f.particles.values())
            graph_recalc |= set(filter(lambda x: x is not None, f.opt.values()))
        # add neighbors since those can be changed too (opt/not atoms)
        graph_recalc = self.populate_neighbors(set(product_particles))

        # [rx_break]
        self.process_break(frags, rx)
        # [rx_update]
        self.process_update(frags, rx)
        # (only non skipped get passed here) remove non graph fragment reactants
        for frag in frags:
            if frag.graph is None:
                self.remove_fragment(frag)

        # graphs get recalculated later over the same particles
        self.remove_overlapping_graphs(graph_recalc)

        frags = Fragments.from_frags(frags)
        # change the interactions - TODO fix it, allow choosing the particles
        for product, selection in rx.products:
            selected = [frags.get((i, name)) for i, name in selection]
            if len(selection) == 0:
                selected = product_particles
            frag = self.type_lookup.get(product)  # frag type
            if frag is None:
                raise ValueError(f"Can't find mol {product}.")
            elif not isinstance(frag, MolFragment):
                raise ValueError(f"Product {product}, is"
                                 " not a mol fragment type."
                                 f" It is: {frag}"
                                 "Use [frag_from] to create fragments "
                                 "in reactions.")

            self.instantiate_over_existing(frag, selected, reaction=True)
        self.instantiate_over_existing(rx.product, frags, template=rx)

        self.try_match_graphs(graph_recalc)

    def post_modification(self, i: int) -> None:
        # hook that only gets called after modification
        for reporter in self.reporters:
            reporter.post_modification(i)

    def detection_modification(self, i: int) -> int:
        self.pre_detection(i)
        init_map = self.get_init_map()
        pos, box = self.system.get_positions()
        query, frag_to_nlist, nlist_to_frag = self.get_neighbor_list(box, pos)
        # tree of frag combinations to check
        reactions, skip, _ = self.detection_over_types(
            [], set(), self.get_reaction_tree(), [], [], init_map, pos, box,
            query, frag_to_nlist, nlist_to_frag
        )

        if len(reactions) == 0:
            return 0
        self.pre_modification(reactions, i)
        for (frags, rx) in reactions:
            self.modification(frags, rx)
        self.post_modification(i)
        return len(reactions)
