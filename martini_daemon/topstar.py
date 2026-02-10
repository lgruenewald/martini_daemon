import random
from .sysstar import SysStar
from .forces.force import Interaction
from .reporters.reporter import Reporter
from .graph import Graph, GraphMatch, GraphAtomType, match_atoms
from .fragment import Fragment
from .reaction_template import ReactionTemplate
from .molecule import Molecule
from .detection import detection_all
#from .detection2 import detection_all
import numpy as np


class TopStar():
    # ======= (1/3) Building things =======
    def __init__(self, system, logger, nlist_cutoff, max_absolute_rate,
                 highest_probability, respos):
        # == 1. Fragments ==
        # TODO rework this into a dense structure
        # fragment by unique id
        self.frag_list: dict[int, Fragment] = {}
        # next unique id
        self.next_frag_id: int = 0
        # TODO rework this into a dense structure
        # for every atom_id have a list of fragments it is in
        # and a list of interactions it is in
        self.defrag_list: list[list[int]] = []
        # TODOish dense structure?
        self.interaction_list: list[list[Interaction]] = []
        # set by self.update_frag_counts()
        # TODO too much state synchronization...
        self.frag_counts: dict[str, int] = {}

        # == 2. Graphs ==
        self.graphs: dict[str, Graph] = {}

        # == 3. Molecules ==
        self.molecules: dict[str, Molecule] = {}

        # == 4. Reactions ==
        # (rx1, rx2, rx3, rx4) -> templates for this reactant combo
        self.reactions: dict[
            tuple[str, str | None, str | None, str | None],
            list[ReactionTemplate]
        ] = {}
        # are there any reactions that continue a certain combo of reactants
        self.r_continue: dict[tuple[str, str | None, str | None], bool] = {}
        # which atom to use within a fragment type for the neighbor listt
        self.neighbor_atom_map: dict[str, int] = {}

        # == 5. Misc ==
        self.reporters: list[Reporter] = []
        self.system: SysStar = system
        self.logger = logger
        self.nlist_cutoff: float = nlist_cutoff
        # warmup phase for rate control
        self.absolute_rate: float = -1.
        self.max_absolute_rate: float = max_absolute_rate
        self.highest_probability: float = highest_probability
        self.respos = respos
        # a list of what was in [molecules] in the .top,
        # used for analysis only at the moment
        self.initial_molecules: list[tuple[str, int, int]] = []
        # set in self.init_dm
        self.out_name = None
        random.seed()

    # Save/load helpers - TODO make load an option in the constructor
    def save(self, f):
        """
            Serializes T* into bytes, writes it to f

            Data saved:
            - frag_list
            - graphs (graph fragment list)
            - molecules (molecule definition list)
            - reactions (reaction template list)
        """
        # relies on pickling to handle references to other classes right
        f.dump(random.getstate())

        f.dump(self.frag_list)
        f.dump(self.graphs)
        f.dump(self.reactions)
        f.dump(self.next_frag_id)
        f.dump(self.defrag_list)
        f.dump(self.r_continue)
        f.dump(self.neighbor_atom_map)
        f.dump(self.absolute_rate)
        f.dump(self.initial_molecules)
        f.dump(self.respos)

        # unpicklable because they reference force
        f.dump(len(self.molecules))
        for k, v in self.molecules.items():
            f.dump(k)
            assert k == v.molecule_name
            v.save(f)

    def load(self, f):
        """
            Deserializes bytes (read from f) into T* (self)

            Should be called after S* is deserialized and set, because
            of the interaction list.
        """
        random.setstate(f.load())

        self.frag_list = f.load()
        self.graphs = f.load()
        self.reactions = f.load()
        self.next_frag_id = f.load()
        self.defrag_list = f.load()
        self.r_continue = f.load()
        self.neighbor_atom_map = f.load()
        self.absolute_rate = f.load()
        self.initial_molecules = f.load()
        self.respos = f.load()

        for i in range(f.load()):
            k = f.load()
            v = Molecule(k)
            v.load(f, self.system)
            self.molecules[k] = v

        # interaction list gets loaded different because we can't pickle
        # openmm things -> can't pickle Force objects
        #
        # but we can assume that everything added to any Force would have
        # resulted in an interaction, except non bonded
        #
        # though this is still fragile code and should be improved later

        self.interaction_list = [[] for _ in range(self.system.len_atoms())]
        for force in self.system.modular_forces:
            if force == self.system.nonbonded_force:
                continue
            for i in range(len(force)):
                inter = Interaction(force, i)
                # TODO make this better
                if force._list[i] is None:
                    continue
                members = inter.get_members()
                for member in members:
                    self.interaction_list[member].append(inter)

    def add_reporter(self, reporter) -> None:
        self.reporters.append(reporter)

    def new_molecule(self, name: str) -> Molecule:
        if self.molecules.get(name):
            raise ValueError(f"Second definition of molecule type {name}")
        molecule = Molecule(name)
        self.molecules[name] = molecule
        return molecule

    def new_graph(self, graph: Graph) -> None:
        self.graphs[graph.name] = graph

    # TODO fix this monster
    def new_reaction(self, reaction: ReactionTemplate) -> Molecule:
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
                    raise ValueError(
                        "r_max with index higher than the number of reactants"
                    )
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
                graph1 = self.graphs[name1]
                graph2 = self.graphs[name2]
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
        return self.molecules[reaction.name]

    # Graph helpers

    def try_match_graphs(self, atoms: set[int], molname=None) -> None:
        """
        Given a set of atoms, find all graph matches of all known
        graphs and add them to the fragment list.

        Should be called after interactions (self.interaction_list) have
        been updated.

        Initial system construction: call this for every molecule, specify
        molname.

        Reaction: Atoms should be all affected atoms in a reaction,
        and all their neighbors. Remove all fragments that contain atoms
        on which try_match_graphs is called first.

        If molname is specified, it means that we can be assured that
        the same matches are going to happen when we call it with the
        same molname again, so we can cache the results and speed up
        this function call later.
        """

        matches: list[GraphMatch] = []
        for graph in self.graphs.values():
            # for each possible graph to match
            # don't match if it's only the specific molecule
            if len(graph.molecules) > 0:
                if molname not in graph.molecules:
                    continue
            matches += match_atoms(
                graph, atoms, self.system,
                self.interaction_list
            )

        for m in matches:
            inst = Fragment(m.graph, self.next_frag_id)
            self.frag_list[self.next_frag_id] = inst
            self.next_frag_id += 1

            for key, _, _, type in m.graph.atoms:
                val = m.atoms.get(key)
                # key - name in the graph
                # val - atom id
                if val is None:
                    inst.atoms.append(-1)
                else:
                    inst.atoms.append(val)
                    self.defrag_list[val].append(inst.frag_id)

    def instantiate(self, molname: str) -> None:
        """
        Adds a fresh new copy of moleculetype molname to the system.
        R
        """
        # molname must refer to a Molecule type
        mol = self.molecules.get(molname)
        if mol is None:
            raise ValueError(f"Can't find mol {molname}")

        # add atoms to S* and other bookkeeping
        atoms = []
        prev_resnum = 0
        for atom in mol.atoms:
            type, resnum, resname, atomname, chargegr, charge, mass = atom
            if resnum != prev_resnum:
                self.system.new_residue()
                prev_resnum = resnum
            p = self.system.add_atom(atomname, resname, type, charge, mass, 1, 0.5)
            atoms.append(p)
            self.defrag_list.append([])
            self.interaction_list.append([])

        # instantiate interactions
        self.instantiate_over_existing(mol, atoms)

        # add graphs to system
        self.try_match_graphs(set(atoms), molname)

        # do initial molecules info, used e.g. in helpers/monomer
        if (
            len(self.initial_molecules) > 0
            and self.initial_molecules[-1][0] == molname
        ):
            _, n, n_atoms = self.initial_molecules[-1]
            assert len(atoms) == n_atoms
            self.initial_molecules[-1] = (molname, n+1, n_atoms)
        else:
            self.initial_molecules.append((molname, 1, len(atoms)))

    def instantiate_over_existing(
        self, mol: Molecule, atoms: list[Fragment] | list[int],
    ) -> None:
        """
        Takes a name of a molecule, adds interactions to those atoms
        according to the molecule.
        """
        def index_pair(index: int | tuple[int, int]) -> int:
            # abstraction, because of the duality of Molecule
            if type(index) is int:
                return atoms[index]
            else:
                idi, atomi = index
                return atoms[idi].atoms[atomi]

        # exclusions
        for (i, j) in mol.exclusions:
            pi = index_pair(i)
            pj = index_pair(j)
            if pi == -1 or pj == -1:
                # optional atom missing
                continue
            if pi < pj:
                e = self.system.exclusions.add(pi, pj)
                self.interaction_list[pi].append(e)
                self.interaction_list[pj].append(e)
        # posres
        for (i, kx, ky, kz) in mol.posres:
            if self.respos is None:
                raise ValueError(
                    "respos is None but there are position restraints."
                )
            pi = index_pair(i)
            x0, y0, z0 = self.respos[pi]
            self.system.posres.add((pi), (kx, ky, kz, x0, y0, z0))
        # generic interactions
        for (force, members, params) in mol.interactions:
            member_atoms = []
            missing_opt = False
            for x in members:
                atom = index_pair(x)
                if atom == -1:
                    missing_opt = True
                    break
                member_atoms.append(atom)
            if missing_opt:
                # optional missing => skip
                continue
            f = force.add(member_atoms, params)
            for member in member_atoms:
                self.interaction_list[member].append(f)

    def update_frag_counts(self):
        """
        Sets self.frag_counts. Call after every modification.

        self.frag_counts is used notably for rate calculation.
        """
        # TODO too much state sync..
        self.frag_counts = {}
        for v in self.frag_list.values():
            k = v.name
            if self.frag_counts.get(k) is None:
                self.frag_counts[k] = 1
            else:
                self.frag_counts[k] += 1

    # ======= (2/3) Detection things - see detection.pyx =======
    def init_dm(self, name) -> None:
        # TODO remove this function, reduce fragility that way
        # called exactly once after parsing or loading from file is finished
        self.update_frag_counts()
        self.out_name = name
        for reporter in self.reporters:
            reporter.init_dm(name)

    def update_observed_rate(self, rx: ReactionTemplate):
        # update_observed_rate only called if it is defined => not none
        # relative_rate is parsed as "positive" => non zero, safe to divide
        rate = rx.reaction_counter / rx.relative_rate
        for reactant in rx.reactants:
            if self.frag_counts.get(reactant) in {0, None}:
                # no reactant => no reaction
                # keeps old observed_rate and we avoid divisions by 0
                return
            rate /= self.frag_counts[reactant]

        # if it's not initialized yet make an initial value
        if rx.observed_rate is None:
            rx.observed_rate = rate
        else:
            # otherwise smooth it
            rx.observed_rate = rate * 0.01 + rx.observed_rate * 0.99

        predicted = rx.observed_rate * self.highest_probability

        # if the prediction value goes below 0, 0 it and issue a warning
        if predicted < 0.:
            self.logger.warn(
                f"The smoothed rate for reaction {rx.name} went below 0."
            )
            rx.observed_rate = 0.
            predicted = 0.

        # set the relative rate = 1 value to the slowest reaction
        if self.absolute_rate is None or predicted < self.absolute_rate:
            self.absolute_rate = predicted

    def detection(
        self, step: int, box, pos
    ) -> list[tuple[list[Fragment], ReactionTemplate]]:
        for reporter in self.reporters:
            reporter.pre_detection(step, self.out_name)
        rxs = []
        for more_rxs in self.reactions.values():
            for rx in more_rxs:
                rxs.append(rx)
        reactions = detection_all(
            self, box, pos
        )
        """
        reactions = detection_all(
            self.frag_list, rxs,
            box, pos,
            self.nlist_cutoff,
            self.absolute_rate,
            # TODO save state (not sure there is much point)
            np.random.default_rng()
        )
        """
        # preparations for the next step
        self.absolute_rate = self.max_absolute_rate
        for _, rxs in self.reactions.items():
            for rx in rxs:
                if rx.relative_rate is not None:
                    self.update_observed_rate(rx)
                rx.reaction_counter = 0
        if self.absolute_rate is None:
            self.absolute_rate = -1.
        return reactions

    # ======= (3/3) Modification things =======
    def clean_defrag(self, frag: Fragment, atom: int) -> None:
        """Removes references to frag in the defrag_list for atom_id atom."""
        self.defrag_list[atom] = list(filter(
            lambda x: x != frag.frag_id,
            self.defrag_list[atom]
        ))

    def remove_fragment(self, frag: Fragment) -> None:
        """Removes a fragment from frag_list and defrag_list"""

        for atom in frag.atoms:
            if atom != -1:
                self.clean_defrag(frag, atom)

        del self.frag_list[frag.frag_id]

    def remove_interaction(self, interaction: Interaction) -> None:
        """Removes an interaction from interaction_list and S*"""
        for atom in interaction.get_members():
            self.interaction_list[atom] = list(filter(
                lambda x: x != interaction,
                self.interaction_list[atom]
            ))
        interaction.remove()

    def process_break(self, frags: list[Fragment], rx: ReactionTemplate,
                      ) -> None:
        """Process [break] in rx over frags."""

        for group in rx.break_groups:
            group_atoms = []
            all_found = True
            for id, atom in group:
                atom = frags[id].atoms[atom]
                if atom == -1:
                    all_found = False
                    break
                group_atoms.append(atom)
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
        """Process [rx_update] in rx over frags."""

        for group in rx.update_groups:
            group_atoms = []
            all_found = True
            for id, atom in group:
                atom = frags[id].atoms[atom]
                if atom == -1:
                    all_found = False
                    break
                group_atoms.append(atom)
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

    def populate_neighbors(self, atoms: set[int]) -> set[int]:
        """
        For a set of atoms, return a set that also contains their neighbors.

        Neighbor = shared interaction (e.g. bond, angle, exclusion...)
        """
        res = set()
        for atom in atoms:
            res.add(atom)
            for inter in self.interaction_list[atom]:
                for member in inter.get_members():
                    res.add(member)
        return res

    def remove_overlapping_graphs(self, atoms: set[int]) -> None:
        """Remove all fragments in T* that contain any of the given atoms."""
        for atom in atoms:
            for frag_id in self.defrag_list[atom][:]:
                frag = self.frag_list.get(frag_id)
                if frag is not None and frag.graph is not None:
                    self.remove_fragment(frag)

    def template_update_atoms(
        self, frags: list[Fragment], rx: ReactionTemplate
    ) -> None:
        """Process [rename], [retype], [recharge], [remass]."""
        # re* overrides everything
        for (id, atom, new_name) in rx.renames:
            atom_id = frags[id].atoms[atom]
            if atom_id == -1:
                continue
            self.system.rename(atom_id, new_name)
        for (id, atom, new_type) in rx.retypes:
            atom_id = frags[id].atoms[atom]
            if atom_id == -1:
                continue
            self.system.retype(atom_id, new_type)
        for (id, atom, new_charge) in rx.recharges:
            atom_id = frags[id].atoms[atom]
            if atom_id == -1:
                continue
            self.system.recharge(atom_id, new_charge)
        for (id, atom, new_mass) in rx.remasses:
            atom_id = frags[id].atoms[atom]
            if atom_id == -1:
                continue
            self.system.remass(atom_id, new_mass)

    def modification(
        self,
        i: int,
        reactions: list[tuple[list[Fragment], ReactionTemplate]]
    ) -> None:
        """Modification helper for the D/M algorithm
        """

        for reporter in self.reporters:
            reporter.pre_modification(i, reactions, self.out_name)

        completed_reactions = []

        for (frags, rx) in reactions:
            # check if all reactants still exist
            if any(
                [self.frag_list.get(frag.frag_id) is None for frag in frags]
            ):
                # we have no way of detecting
                # this in the detection algorithm, since this is
                # about graphs being recalculated during reactions
                # potentially invalidating other graphs that react
                # in the same frame
                continue

            # just normal atoms to instantiate products over
            product_atoms = []
            # which atoms to recalculate graphs over
            graph_recalc = set()
            for f in frags:
                product_atoms += list(filter(lambda x: x != -1, f.atoms))
                graph_recalc |= set(filter(lambda x: x != -1, f.atoms))
            # add neighbors since those can be changed too (opt/not atoms)
            graph_recalc = self.populate_neighbors(set(product_atoms))

            # [rx_break]
            self.process_break(frags, rx)
            # [rx_update]
            self.process_update(frags, rx)

            # graphs get recalculated later over the same atoms
            self.remove_overlapping_graphs(graph_recalc)

            self.template_update_atoms(frags, rx)
            self.instantiate_over_existing(self.molecules[rx.name], frags)

            self.try_match_graphs(graph_recalc)

            completed_reactions.append((frags, rx))

        self.update_frag_counts()

        for reporter in self.reporters:
            reporter.post_modification(i, self.out_name, completed_reactions)

        return completed_reactions

    def toggle_sc(
            self, reactions: list[tuple[list[Fragment], ReactionTemplate]],
            toggle: bool
    ) -> None:
        for (frags, rx) in reactions:
            for (id, atom, new_lam, new_alpha) in rx.soft_core:
                atom_id = frags[id].atoms[atom]
                if atom_id == -1:
                    continue
                if toggle:
                    self.system.update_sc(atom_id, new_lam, new_alpha)
                else:
                    self.system.update_sc(atom_id, 1, 0.5)

#            for f in frags:
#                for atom_id in f.atoms:
#                    if atom_id == -1:
#                         continue
#                    if toggle:
#                        self.system.update_sc(atom_id, 0.6, 0.5)
#                    else:
#                        self.system.update_sc(atom_id, 1, 0.5)
        return

