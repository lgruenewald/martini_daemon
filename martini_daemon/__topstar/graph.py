# helper function for T* graph based fragments
from __future__ import annotations

from enum import Enum
from fnmatch import fnmatch

from ..__core import BondedForce, System
from ..__parser import ParseException


# === MAIN CLASSES ===
class GraphAtomType(Enum):
    NORMAL = 0
    OPT = 1
    NOT = 2


class Graph:
    def __init__(self, name: str | None) -> None:
        """
        Create a Graph instance.

        The class constructed from [graph]/[frag] directives that contains all
        the information the user provided about a graph.

        Note: Graphs contain atoms that each have a name (which is what's used
        in input files) and an index (used only internally, may show up in some outputs).
        The index is used for canonical ordering in fragments, and name->index
        is also resolved during parsing in modification template's
        parse_pair and parse_index.

        :param name: Name of the graph.

        Attributes:
        - name - graph name
        - atoms - list of (graph_atom_name, name_pattern, type_patter, graph_atom_type)
        - interactions - list of (interaction_type, list of graph_atom_name)
        - equivalents - list of sets of str
        - atom_name_to_index - dict of str and int

        """
        self.name: str | None = name
        self.atoms: list[tuple[str, str, str, GraphAtomType]] = []
        self.interactions: list[tuple[str, list[str]]] = []
        self.equivalents: list[set[str]] = []
        self.atom_name_to_index: dict[str, int] = {}

    def finish_init(self) -> None:
        """
        Must be called after parsing the graph and before it's used.

        Validates graphs and errors on malformed graphs.
        Raises ParseExceptions.
        """
        if self.name is None:
            raise ParseException("Graph has no name.")
        nodes: dict[str, set[str]] = {}
        if len([x for x in self.atoms if x[3] == GraphAtomType.NORMAL]) == 0:
            raise ValueError("Graph must contain at least one normal atom")
        if self.atoms[0][3] != GraphAtomType.NORMAL:
            raise ParseException(
                "First atom in graph must be a normal atom (not opt or not)."
            )
        for i, (name, _, _, _) in enumerate(self.atoms):
            if nodes.get(name) is not None:
                raise ParseException(
                    f"Same graph has multiple atoms of the same name: {name}."
                )
            nodes[name] = set()
            self.atom_name_to_index[name] = i
        for filter_str, atoms in self.interactions:
            if len(atoms) < 2:
                raise ParseException(
                    f"Interaction {filter_str} must have at least "
                    f"two atom members. Only found {atoms}."
                )
            num_special = 0
            for atom in atoms:
                if self.atoms[self.atom_name_to_index[atom]][3] != GraphAtomType.NORMAL:
                    num_special += 1
                if nodes.get(atom) is None:
                    raise ParseException(
                        f"Interaction {filter_str} for atoms {atoms} "
                        f"references undefined atom name {atom}."
                    )
                for other_atom in atoms:
                    if atom != other_atom:
                        nodes[atom].add(other_atom)
            if num_special > 1:
                raise ParseException(
                    f"Interaction {filter_str} for atoms {atoms} "
                    + "references more than one optional or forbidden atoms. "
                    + "Grouping forbidden/optional atoms is not supported."
                )

        # if there is more than 1 atom, they all must be connected to at least
        # one other atom with an interaction
        marked: set[str] = set()

        # go from atom 1 and mark all
        def mark(atom_: str) -> None:
            if atom_ in marked:
                return
            marked.add(atom_)
            for other_atom in nodes[atom_]:
                mark(other_atom)

        mark(self.atoms[0][0])
        if len(marked) != len(self.atoms):
            assert len(self.atoms) == len(nodes.keys())
            missing = set(nodes.keys()) - marked
            raise ParseException(
                f"Invalid graph. The graph is not all connected to itself. Add"
                f" interactions to connect it. Disconnected atom: {missing}."
            )


# === MATCH HELPERS ===
class GraphMatch:
    # mapping of atoms -> atom_id
    graph: Graph
    atoms: dict[str, int]
    # interaction matches, lists of None or Interaction
    interactions: list[None | tuple[str, int]]
    # reverse of atoms
    rev_atoms: dict[int, str]
    matched_inter: set[tuple[str, int]]
    next_inter: int

    def __init__(self, graph: Graph) -> None:
        """
        Create a Helper class that represents a (partially) mapped out graph to S*.

        Attributes:
        - graph - reference to a Graph instance that was (partially) matched
        - atoms - dict of atoms that were matched, dict of atom name to atom_id
        - interactions - interaction matches for each interaction in graph. All are None or an Interaction
            (a tuple of force name and bond_id).
        - rev_atoms - atom_id to atom name mapping
        - matched_inter - same as interactions, but as a set and excluding None's
        - next_inter - internal state for the graph matching algorithm, it should point to the next
            interaction index in the graph that was not attempted to be filled yet.

        """
        self.graph = graph
        self.atoms = {}
        self.interactions = [None for _ in graph.interactions]
        self.rev_atoms = {}
        self.matched_inter = set()
        self.next_inter = 0

    def copy(self) -> GraphMatch:
        res = GraphMatch(self.graph)
        res.atoms = self.atoms.copy()
        res.interactions = self.interactions.copy()
        res.rev_atoms = self.rev_atoms.copy()
        res.matched_inter = self.matched_inter.copy()
        # next_inter intentionally not copied
        return res

    def add_atom(self, name: str, atom_num: int) -> None:
        """Add an graph atom name <=> atom_id mapping to the partial match."""
        assert self.atoms.get(name) is None and atom_num not in self.rev_atoms
        self.atoms[name] = atom_num
        self.rev_atoms[atom_num] = name

    def add_inter(self, id: int, inter: tuple[str, int]) -> None:
        assert self.interactions[id] is None and inter not in self.matched_inter
        self.interactions[id] = inter
        self.matched_inter.add(inter)

    def is_complete(self) -> bool:
        # interactions are checked during construction
        # only partials that this should be called on are ones where
        # every time a new atom is added, all interactions are enforced
        # that can be enforced when adding that atom
        # TLDR only checks atoms, not interactions
        for name, _, _, atom_type in self.graph.atoms:
            found = self.atoms.get(name) is not None
            if atom_type == GraphAtomType.NORMAL and not found:
                return False
        return True

    def is_acceptable(self) -> bool:
        # is_complete + NOT type checking
        for name, _, _, atom_type in self.graph.atoms:
            found = self.atoms.get(name) is not None
            if atom_type == GraphAtomType.NORMAL and not found:
                return False
            if atom_type == GraphAtomType.NOT and found:
                return False
        return True

    def is_equal(self, other: GraphMatch) -> bool:
        # note1: also see note for is_complete / is_acceptable -> only checks atoms
        # note2: equivalent atoms in graph are exchangeable
        if self.graph != other.graph:
            return False
        if len(self.atoms) != len(other.atoms):
            return False
        for name, id in self.atoms.items():
            other_name = other.rev_atoms.get(id)
            if other_name is None:
                return False
            if other_name != name and all(
                name not in eqs or other_name not in eqs
                for eqs in self.graph.equivalents
            ):
                return False
        return True


class AtomCache:
    def __init__(self, system: System) -> None:
        """Create a Helper class that groups S* information and provides helper query functions to it."""
        self.system = system

    def neighbors(self, atom: int) -> set[int]:
        res = set()
        inters = self.system.get_interactions_for_atom(atom)
        if len(inters) == 0:
            return res
        for inter in inters:
            for member in self.system.get_members(*inter):
                res.add(member)
        res.remove(atom)
        return res

    def check_atom_interactions(
        self, g_atom: str, atom_id: int, partial: GraphMatch
    ) -> tuple[bool, list[tuple[int, tuple[str, int]]]]:
        """
        Check if atom_id can be g_atom in the graph, based on interaction filters.

        Returns True if adding g_atom=atom_id to the graph match is
        possible (all interaction requirements fulfilled).
        Returns False if there is an interaction requirement violated
        (missing interaction that should be there).

        :param g_atom: Graph atom name
        :param atom_id: S* atom ID candidate for g_name
        :param partial: Partial graph match that g_atom=atom_id is considered
            for.

        Note: assumes g_atom=atom_id is not a part of partial yet.

        Note2: this function does not mutate partial.
        """
        # g_ prefix -> graph things
        # s_ prefix -> S* things
        s_inters = self.system.get_interactions_for_atom(
            atom_id
        )  # interactions for atom_id in S*
        # interactions already considered
        skip: set[tuple[str, int]] = {i for i in partial.interactions if i is not None}
        matches: list[tuple[int, tuple[str, int]]] = []  # new matches
        # iterate over all interaction requirements in graph
        for i, (g_type, g_atoms) in enumerate(partial.graph.interactions):
            # g_type -> the filter text of this interaction (e.g. "connection")
            # g_atoms -> list of graph atom names in the interaction

            # already has an interaction for this one, skip
            if partial.interactions[i] is not None:
                continue
            # atom name not in this interaction, skip
            if g_atom not in g_atoms:
                continue
            # members already there in the partial graph match for this graph
            # inter
            g_members = {
                # g_atom: atom_id, otherwise query partial
                partial.atoms.get(name) if name != g_atom else atom_id
                for name in g_atoms
            }
            # missing atom for this interaction? for now ignore, will be
            # matched once it gets filled
            if None in g_members:
                continue

            # find which interaction in S* corresponds to this to-be-filled
            # by now interaction in the graph
            # Note: unsound code below
            # we eagerly take the first interaction that matches, this might
            # not be what we want (e.g. vsite 1 2 3 and vsite 1 2, we might
            # accept vsite 1 2 3 for a constraint vsite 1 2)
            #
            # but this is rarely a problem and would add a lot of complexity
            # for this edge case
            found = False
            for s_inter in s_inters:
                # S* interaction does not fulfill graph type filter
                # we assume that it's a BondedForce
                force = self.system.get_force(s_inter[0])
                assert isinstance(force, BondedForce)
                if not force.passes_filter(g_type):
                    continue
                # s_inter already used
                if s_inter in skip:
                    continue
                # S* Interaction members
                s_members = set(self.system.get_members(*s_inter))
                if len(g_members - s_members) > 0:
                    # S* can contain extra members, but all Graph ones
                    # should be fulfilled
                    continue
                # got here? matching S* inter = Graph inter
                found = True
                matches.append((i, s_inter))
                skip.add(s_inter)
                break
            if not found:
                return False, []
        return True, matches

    def is_name_type(self, name_filter: str, type_filter: str, atom_id: int) -> bool:
        name = self.system.get_name(atom_id)
        type_ = self.system.get_type(atom_id)
        return fnmatch(name, name_filter) and fnmatch(type_, type_filter)


# === MAIN MATCHING ALGO ===
def match_atoms(graph: Graph, atoms: set[int], system: System) -> list[GraphMatch]:
    """Return all unique graph matches for graph against a given set of atoms."""
    # 1. build a atom cache
    cache = AtomCache(system)
    queue: list[GraphMatch] = []  # partial matches
    results: list[GraphMatch] = []  # complete matches
    # 2. find starting matches of a single atom
    # find all normal/optional single atom matches and add them to the queue
    for atom in atoms:
        for name, name_filter, type_filter, graph_type in graph.atoms:
            if graph_type is GraphAtomType.NOT:
                continue
            if cache.is_name_type(name_filter, type_filter, atom):
                new_match = GraphMatch(graph)
                new_match.add_atom(name, atom)
                queue.append(new_match)

    # main / queue loop, formerly it was using recursion
    while len(queue) > 0:
        cmatch = queue.pop(0)
        # 1. find one interaction that is missing but already has at least one
        # atom
        assert len(cmatch.interactions) == len(graph.interactions)
        any_found = False

        # go through the remaining missing interactions in this graph match
        # the first suitable is picked, and all atom candidates for it
        # get added to the queue with next_inter reset.
        while cmatch.next_inter < len(cmatch.interactions):
            inter_id = cmatch.next_inter
            cmatch.next_inter += 1
            inter = cmatch.interactions[inter_id]
            if inter is not None:
                # already filled, skip
                continue
            # corresponds to a inter in the graph (e.g. connection cc1 cc2)
            inter_filter, inter_atoms = graph.interactions[inter_id]
            # find one of the atoms in the interaction
            missing: list[str] = []  # <- holes in the graph (graph atom names)
            len_filled = 0  # <- # of filled atoms
            last_filled: int  # <- one of the atoms already in the inter
            for atom in inter_atoms:
                atom_id = cmatch.atoms.get(atom)
                if atom_id is None:
                    missing.append(atom)
                else:
                    len_filled += 1
                    last_filled = atom_id
            # no atom yet, so we can't begin to fill it, skip for now
            if len_filled == 0:
                continue
            # 2. query ALL possible new matches in the interaction and put
            # them in the queue
            # new matches have to:
            # a) match name/type
            # b) pass check_atom_interactions
            neighbors = cache.neighbors(last_filled)
            for new_atom in neighbors:
                # only look for new atoms
                if cmatch.rev_atoms.get(new_atom) is not None:
                    continue
                # exhaustively treat all matches here: all holes with all
                # neighbors
                for cmissing in missing:
                    cmissing_id = graph.atom_name_to_index[cmissing]
                    _, name_filter, type_filter, _ = graph.atoms[cmissing_id]
                    if not cache.is_name_type(name_filter, type_filter, new_atom):
                        continue
                    valid, matches = cache.check_atom_interactions(
                        cmissing, new_atom, cmatch
                    )
                    if not valid:
                        continue
                    new_cmatch = cmatch.copy()
                    new_cmatch.add_atom(cmissing, new_atom)
                    for i, inter in matches:
                        new_cmatch.add_inter(i, inter)
                    queue.append(new_cmatch)
                    # we are exhaustive with filling just this one gap in the
                    # graph so we can commit to the set of things in the queue
                    # added here, and we no longer need cmatch after this
                    # inner loop
                    any_found = True

            # 3. break because we only try to progress on one
            # interaction per queue loop
            break

        # assuming next_inter still points to an interaction,
        # if any_found is False, readd to the queue (next_inter has been incr)
        # this is how we query all interactions one by one, even if we can't
        # fill it all immediately
        if not any_found:
            # cmatch was not consumed
            if cmatch.next_inter < len(cmatch.interactions):
                # there is more to check
                queue.append(cmatch)
            else:
                # if next_inter is already the last one, check if it is
                # complete, acceptable and non duplicate.
                # If all three add to result list.
                if cmatch.is_complete() and cmatch.is_acceptable():
                    duplicate = False
                    for res in results:
                        if cmatch.is_equal(res):
                            duplicate = True
                            break
                    if not duplicate:
                        results.append(cmatch)

    return results
