import random
from .sysstar import SysStar
from .forces.force import Interaction
from .reporters.reporter import Reporter
from .graph import GraphFragment, GraphMatch, GraphAtomType, match_particles
from .fragment import Fragment
from .reaction_template import ReactionTemplate
from .mol_fragment import MolFragment
from .detection import detection_all

random.seed()

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
    graph_fragment_map: dict[str, GraphFragment]

    # type name -> type, used for instantiation
    # TODO only use it for MolFragment
    type_lookup: dict[str, MolFragment]

    # reactant type -> reactions map
    reactions: dict[tuple[str, str, str, str], list[ReactionTemplate]]
    # reactant type -> are there any higher order reactions with this combo map
    r_continue: dict[tuple[str, str, str], bool]
    # reactant type -> neighbor list atom index map
    neighbor_atom_map: dict[str, int]

    system: SysStar

    reporters: list[Reporter]

    def __init__(self, system, logger, nlist_cutoff):
        # TODO better suite the bookkeeping to the algorithm used and remove
        # the other ones
        self.frag_list = {}
        self.next_frag_id = 0
        self.defrag_list = []
        self.interaction_list = []
        self.frag_fragments = []
        self.mol_fragments = []
        self.reactions: dict[tuple[str, str, str, str], list[ReactionTemplate]] = {}
        self.reaction_list: list[ReactionTemplate] = []
        self.r_continue = {}
        self.neighbor_atom_map = {}
        self.type_lookup = {}
        self.subfrag_map = {}
        self.reporters = []
        self.graph_fragment_list = []
        self.graph_fragment_map = {}
        self.system = system
        self.logger = logger
        self.nlist_cutoff = nlist_cutoff

    def add_reporter(self, reporter) -> None:
        self.reporters.append(reporter)

    def new_mol_fragment(self, name: str) -> MolFragment:
        if self.type_lookup.get(name):
            raise ValueError(f"Second definition of fragment type {name}")
        mol_fragment = MolFragment(name)
        self.type_lookup[name] = mol_fragment
        return mol_fragment

    def new_graph_fragment(self, frag: GraphFragment) -> None:
        # TODO unify
        self.graph_fragment_list.append(frag)
        self.graph_fragment_map[frag.name] = frag

    def new_reaction(self, reaction: ReactionTemplate) -> MolFragment:
        # reactant names
        r1 = reaction.reactants[0]
        r2 = reaction.reactants[1] if len(reaction.reactants) >= 2 else None
        r3 = reaction.reactants[2] if len(reaction.reactants) >= 3 else None
        r4 = reaction.reactants[3] if len(reaction.reactants) >= 4 else None
        key = (r1, r2, r3, r4)
        # r_max validation
        if len(reaction.reactants) >= 2:
            r_max_connected = []
            for i in range(len(reaction.reactants)):
                r_max_connected.append({i})
            for i1, aindex1, i2, aindex2, dist in reaction.distance_max:
                if i1 > len(reaction.reactants) or i2 > len(reaction.reactants):
                    raise ValueError("r_max with index higher than the number of reactants")
                name1 = key[i1]
                name2 = key[i2]
                if i1 == i2:
                    continue
                # distance warning
                if dist * 2. > self.nlist_cutoff:
                    self.logger.warn(
                        f"Warning: r_max for reaction {reaction.name}"
                        f" has a r_max between reactant {i1} and {i2}"
                        f" atoms {aindex1} {aindex2} of {dist},"
                        f" which is too large relative to the"
                        f" neighborlist cutoff of {self.nlist_cutoff}."
                    )
                # graph 1 and 2
                graph1 = self.graph_fragment_map[name1]
                graph2 = self.graph_fragment_map[name2]
                # index to type
                _, _, _, type1 = graph1.atoms[aindex1]
                _, _, _, type2 = graph2.atoms[aindex2]
                # we only care about r_max that's guaranteed to be there
                if not (type1 == type2 and type1 == GraphAtomType.NORMAL):
                    continue
                # save this as a valid r_max
                r_max_connected[i1].add(i2)
                r_max_connected[i2].add(i1)
                _, dist1 = self.neighbor_atom_map.get(name1) or (0, 0.)
                _, dist2 = self.neighbor_atom_map.get(name2) or (0, 0.)
                if dist > dist1:
                    self.neighbor_atom_map[name1] = (aindex1, dist)
                if dist > dist2:
                    self.neighbor_atom_map[name2] = (aindex2, dist)
            # graph traversal to see all reactants have a distance cutoff
            marked = set()

            def mark(n):
                if n in marked:
                    return
                marked.add(n)
                for i in r_max_connected[n]:
                    mark(i)
            mark(0)
            if len(marked) != len(reaction.reactants):
                marked_p1 = {x+1 for x in marked}
                raise ValueError("Not all reactants are connected via an r_max"
                                 f"condition of non-optional atoms, "
                                 " therefore reaction "
                                 f"{reaction.name} is invalid. "
                                 f"{marked_p1} are connected to reactant 1.")
        # reactant info buildup
        if r2 is not None:
            self.r_continue[(r1, None, None)] = True
        if r3 is not None:
            self.r_continue[(r1, r2, None)] = True
        if r4 is not None:
            self.r_continue[(r1, r2, r3)] = True

        if self.reactions.get(key) is None:
            self.reactions[key] = []
        self.reactions[key].append(reaction)
        self.reaction_list.append(reaction)
        # reaction product buildup
        if self.type_lookup.get(reaction.name):
            raise ValueError(
                f"Second definition of molfragment type {reaction.name}"
            )
        self.type_lookup[reaction.name] = reaction.product
        return reaction.product

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
            matches += match_particles(
                graph, particles, self.system,
                self.interaction_list
            )

        for m in matches:
            inst = Fragment(m.graph, self.next_frag_id)
            self.frag_list[self.next_frag_id] = inst
            self.next_frag_id += 1

            for key, _, _, type in m.graph.atoms:
                val = m.atoms.get(key)
                # key - name in the graph
                # val - particle id
                if val is None:
                    inst.atoms.append(-1)
                else:
                    inst.atoms.append(val)
                    self.defrag_list[val].append(inst.frag_id)

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
        res = self.instantiate_over_existing(frag, parts)
        # add graphs to system
        self.try_match_graphs(set(parts), frag_name)
        return res

    def instantiate_over_existing(self, molfrag: MolFragment,
                                  frags: list[Fragment] | list[int],
                                  ):
        """Takes a name of a mol fragment, adds interactions to those particles
        according to the mol fragment, or optionally a reaction template.
        """
        def index_pair(frags: list[int] | list[Fragment], index: int | tuple[int, int]) -> int:
            if type(index) is int:
                return frags[index]
            else:
                idi, atomi = index
                return frags[idi].atoms[atomi]
        # TODO good index_pair solution for this one
        # exclusions
        for (i, j) in molfrag.exclusions:
            pi = index_pair(frags, i)
            pj = index_pair(frags, j)
            if pi == -1 or pj == -1:
                # optional atom missing
                continue
            if pi < pj:
                e = self.system.exclusions.add(pi, pj)
                self.interaction_list[pi].append(e)
                self.interaction_list[pj].append(e)
        # generic interactions
        for (force, members, params) in molfrag.interactions:
            member_parts = []
            missing_opt = False
            for x in members:
                part = index_pair(frags, x)
                if part == -1:
                    missing_opt = True
                    break
                member_parts.append(part)
            if missing_opt:
                # optional missing => skip
                continue
            f = force.add(member_parts, params)
            for member in member_parts:
                self.interaction_list[member].append(f)

    # ======= (2/3) Detection things =======
    # see detection.pyx
    def pre_detection(self, step: int) -> None:
        for _, rxs in self.reactions.items():
            for rx in rxs:
                rx.global_counter = 0
        for reporter in self.reporters:
            reporter.pre_detection(step)
    
    def detection(self, step: int, box, pos) -> list[tuple[list[Fragment], ReactionTemplate]]:
        self.pre_detection(step)
        return detection_all(self, box, pos)

    # ======= (3/3) Modification things =======
    def clean_defrag(self, frag: Fragment, part: int) -> None:
        self.defrag_list[part] = list(filter(
            lambda x: x != frag.frag_id,
            self.defrag_list[part]
        ))

    def remove_fragment(self, frag: Fragment) -> None:
        """Removes a fragment from frag_list and defrag_list
        """

        for part in frag.atoms:
            if part != -1:
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
            for id, atom in group:
                part = frags[id].atoms[atom]
                if part == -1:
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
            for id, atom in group:
                part = frags[id].atoms[atom]
                if part == -1:
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

    def template_update_atoms(
                              self, frags: list[Fragment], rx: ReactionTemplate
                             ) -> None:
        # re* overrides everything
        for (id, atom, new_name) in rx.renames:
            part_id = frags[id].atoms[atom]
            if part_id == -1:
                continue
            self.system.rename(part_id, new_name)
        for (id, atom, new_type) in rx.retypes:
            part_id = frags[id].atoms[atom]
            if part_id == -1:
                continue
            self.system.retype(part_id, new_type)
        for (id, atom, new_charge) in rx.recharges:
            part_id = frags[id].atoms[atom]
            if part_id == -1:
                continue
            self.system.recharge(part_id, new_charge)
        for (id, atom, new_mass) in rx.remasses:
            part_id = frags[id].atoms[atom]
            if part_id == -1:
                continue
            self.system.remass(part_id, new_mass)

    def pre_modification(self, rx_list: list[(list, ReactionTemplate)], i
                         ) -> None:
        # hook that gets called after detection, before modification
        # only called if there is any modification going on
        for reporter in self.reporters:
            reporter.pre_modification(rx_list, i)

    def modification(
        self, 
        i: int,
        reactions: list[tuple[list[Fragment], ReactionTemplate]]
    ) -> None:
        """Modification helper for the D/M algorithm
        """

        self.pre_modification(reactions, i)

        for (frags, rx) in reactions:
            # just normal particles to instantiate products over
            product_particles = []
            # which particles to recalculate graphs over
            graph_recalc = set()
            for f in frags:
                product_particles += list(filter(lambda x: x != -1, f.atoms))
                graph_recalc |= set(filter(lambda x: x != -1, f.atoms))
            # add neighbors since those can be changed too (opt/not atoms)
            graph_recalc = self.populate_neighbors(set(product_particles))

            # [rx_break]
            self.process_break(frags, rx)
            # [rx_update]
            self.process_update(frags, rx)

            # graphs get recalculated later over the same particles
            self.remove_overlapping_graphs(graph_recalc)

            self.template_update_atoms(frags, rx)
            self.instantiate_over_existing(rx.product, frags)

            self.try_match_graphs(graph_recalc)

            
        self.post_modification(i)


    def post_modification(self, i: int) -> None:
        # hook that only gets called after modification
        for reporter in self.reporters:
            reporter.post_modification(i)