# helper function for T* graph based fragments
from __future__ import annotations
from enum import Enum
from .forces.force import Interaction
from .sysstar import SysStar
from fnmatch import fnmatch


# === MAIN CLASSES ===
class GraphAtomType(Enum):
    NORMAL = 0
    OPT = 1
    NOT = 2


class Graph():
    """
        The class constructed from [graph]/[frag] directives that contains all
        the information the user provided about a graph.
    """
    name: str
    # only summon this graph on atoms of this molecule during its instantion
    # and only during the start of the simulation (TODO make this more flexible)
    # if empty, all molecules and during reactions too
    molecules: list[str]
    # list[(atom_id, name_pat, type_pat, type)]
    atoms: list[tuple[str, str, str, GraphAtomType]]
    # list[(interaction_type, list[atom_id])]
    interactions: list[tuple[str, list[str]]]
    equivalents: list[set[str]]

    atom_name_to_index: dict[str, int]

    def __init__(self, name):
        self.name = name
        self.molecules = []
        self.atoms = []
        self.interactions = []
        self.equivalents = []
        self.atom_name_to_index = {}

    def finish_init(self):
        """
            Must be called after parsing the graph and before it's used.

            Validates graphs and errors on malformed graphs.
        """
        nodes: dict[str, set[str]] = {}
        if len(list(filter(lambda x: x[3] == GraphAtomType.NORMAL, self.atoms))) == 0:
            raise ValueError("Graph must contain at least one normal atom")
        for i, (name, _, _, _) in enumerate(self.atoms):
            if nodes.get(name) is not None:
                raise ValueError(f"Same graph has multiple atoms of the same name: {name}.")
            nodes[name] = set()
            self.atom_name_to_index[name] = i
        for filter_str, atoms in self.interactions:
            if len(atoms) < 2:
                raise ValueError(f"Interaction {filter_str} must have at least two atom members. Only found {atoms}.")
            for atom in atoms:
                if nodes.get(atom) is None:
                    raise ValueError(f"Interaction {filter_str} for atoms {atoms} references undefined atom name {atom}.")
                for other_atom in atoms:
                    if atom != other_atom:
                        nodes[atom].add(other_atom)

        # if there is more than 1 atom, they all must be connected to at least
        # one other atom with an interaction
        marked: set[str] = set()

        # go from atom 1 and mark all
        def mark(atom: str):
            if atom in marked:
                return
            marked.add(atom)
            for other_atom in nodes[atom]:
                mark(other_atom)
        mark(self.atoms[0][0])
        if len(marked) != len(self.atoms):
            assert len(self.atoms) == len(nodes.keys())
            missing = set(nodes.keys()) - marked
            raise ValueError(f"Invalid graph. The graph is not all connected to itself. Add interactions to connect it. Disconnected atom: {missing}.")


# === MATCH HELPERS ===
class GraphMatch():
    """
        Helper class that represents a (partially) mapped out graph to S*.
    """
    # mapping of atoms -> atom_id
    graph: Graph
    atoms: dict[str, int]
    # interaction matches, lists of None or Interaction
    interactions: list[None | Interaction]
    # reverse of atoms
    rev_atoms: dict[int, str]
    matched_inter: set[Interaction]
    next_inter: int

    def __init__(self, graph: Graph):
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

    def add_atom(self, name: str, atom_num: int):
        assert self.atoms.get(name) is None and atom_num not in self.rev_atoms
        self.atoms[name] = atom_num
        self.rev_atoms[atom_num] = name

    def add_inter(self, id: int, inter: Interaction):
        assert self.interactions[id] is None and inter not in self.matched_inter
        self.interactions[id] = inter
        self.matched_inter.add(inter)

    def is_complete(self):
        # interactions are checked during construction
        # only partials that this should be called on are ones where
        # every time a new atom is added, all interactions are enforced
        # that can be enforced when adding that atom // TLDR only checks atoms
        for (name, _, _, atom_type) in self.graph.atoms:
            found = self.atoms.get(name) is not None
            if atom_type == GraphAtomType.NORMAL and not found:
                return False
        return True

    def is_acceptable(self):
        # is complete + not type checking
        for (name, _, _, atom_type) in self.graph.atoms:
            found = self.atoms.get(name) is not None
            if atom_type == GraphAtomType.NORMAL and not found:
                return False
            if atom_type == GraphAtomType.NOT and found:
                return False
        return True

    def is_equal(self, other: GraphMatch) -> bool:
        # also see note for is_complete / is_acceptable -> only checks atoms
        # equivalent atoms in graph are exchangeable
        if self.graph != other.graph:
            return False
        if len(self.atoms) != len(other.atoms):
            return False
        for name, id in self.atoms.items():
            other_name = other.rev_atoms.get(id)
            if other_name is None:
                return False
            if other_name != name:
                if all(
                    map(
                        lambda eqs: name not in eqs or other_name not in eqs,
                        self.graph.equivalents
                    )
                ):
                    return False
        return True


class AtomCache():
    """
        Helper class that groups S* information and provides helper query
        functions to it.
    """

    def __init__(self, sysstar: SysStar, interactions: list[list[Interaction]]):
        self.sysstar = sysstar
        self.interactions = interactions

    def neighbors(self, atom: int) -> set[int]:
        res = set()
        if len(self.interactions[atom]) == 0:
            return res
        for inter in self.interactions[atom]:
            for member in inter.get_members():
                res.add(member)
        res.remove(atom)
        return res

    def check_atom_interactions(self, atom: str, atom_id: int, partial: GraphMatch) -> bool:
        """
            When adding a new atom, check all interactions that this new atom
            has are complete or still possible.

            Note: partial must not have the atom to be added in it
            yet. atom is the name of the would be added atom, atom_id is the
            index in S*

            Does not mutate partial.
        """
        inters = self.interactions[atom_id]
        skip: set[Interaction] = set()  # interactions already considered
        matches: list[tuple[int, Interaction]] = []
        for i, (inter_type, inter_atoms) in enumerate(partial.graph.interactions):
            # already mapped out, so already full
            if partial.interactions[i] is not None:
                continue
            # atom name not in this interaction, skip
            if atom not in inter_atoms:
                continue
            # members already there in the partial graph match
            inter_g_members = {
                partial.atoms.get(name) if name != atom else atom_id
                for name in inter_atoms
            }
            # missing atom? we don't say anything yet
            if None in inter_g_members:
                continue
            found = False
            for inter in inters:
                # wrong type
                if not inter.is_instance(inter_type):
                    continue
                # already used up
                if inter in skip:
                    continue
                # S* Interaction members
                inter_s_members = set(inter.get_members())
                if len(inter_g_members - inter_s_members) > 0:
                    # S* can contain extra members, but all in graph
                    # should be ones in the graph
                    continue
                # got here? match
                found = True
                matches.append((i, inter))
                skip.add(inter)
                break
            if not found:
                return False, []
        return True, matches

    def is_name_type(self, name_filter: str, type_filter: str, atom_id: int) -> bool:
        name = self.sysstar.get_atom_name(atom_id)
        type, _, _ = self.sysstar.get_atom_details(atom_id)
        return fnmatch(name, name_filter) and fnmatch(type, type_filter)


# === MAIN MATCHING ALGO ===
def match_atoms(
    graph: Graph, atoms: set[int],
    sysstar: SysStar, interactions: list[list[Interaction]]
):
    """
    Return all unique graph matches for graph against a given set of atoms.
    """
    # 1. build a atom cache
    cache = AtomCache(sysstar, interactions)
    queue: list[GraphMatch] = []
    results: list[GraphMatch] = []
    # 2. find starting matches of a single atom
    # any of atoms can be a normal atom obviously
    # optionals included, because it's possible only an optional of the graph
    # is/was in the set of atoms
    # put all of these in a queue
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
        any_found = False
        assert len(cmatch.interactions) == len(graph.interactions)
        for inter_id in range(cmatch.next_inter, len(cmatch.interactions)):
            cmatch.next_inter = inter_id + 1
            inter = cmatch.interactions[inter_id]
            if inter is not None:
                continue
            inter_filter, inter_atoms = graph.interactions[inter_id]
            # find one of the atoms in the interaction
            missing: list[str] = []  # <- holes in the graph (graph atom names)
            len_filled = 0
            last_filled: int  # <- one of the atoms already in the inter
            for atom in inter_atoms:
                atom_id = cmatch.atoms.get(atom)
                if atom_id is None:
                    missing.append(atom)
                else:
                    len_filled += 1
                    last_filled = atom_id
            if len_filled == 0:
                continue
            # by reaching this point we commit to this interaction in this
            # queue iter, increase next_inter
            # 2. query all possible new matches in the interaction and put
            # them in the queue
            # new matches have to:
            # a) match name/type
            # b) pass check_atom_interactions
            neighbors = cache.neighbors(last_filled)
            for new_atom in neighbors:
                # only look for new atoms
                if cmatch.rev_atoms.get(new_atom) is not None:
                    continue
                # exhaustively treat all matches here: all holes with all neighbors
                for cmissing in missing:
                    cmissing_id = graph.atom_name_to_index[cmissing]
                    _, name_filter, type_filter, _ = graph.atoms[cmissing_id]
                    if not cache.is_name_type(name_filter, type_filter, new_atom):
                        continue
                    valid, matches = cache.check_atom_interactions(cmissing, new_atom, cmatch)
                    if not valid:
                        continue
                    new_cmatch = cmatch.copy()
                    new_cmatch.add_atom(cmissing, new_atom)
                    for i, inter in matches:
                        new_cmatch.add_inter(i, inter)
                    queue.append(new_cmatch)
                    # we are exhaustive with filling just this one gap in the graph
                    # so we can commit to the set of things in the queue added here,
                    # and we no longer need cmatch after this inner loop
                    any_found = True

            # note: if any_found is False here, it means there is a hole
            # in the graph that we now know cannot be filled starting at
            # the previous value of cmatch.next_inter

            # 3. break because we only try to progress on one
            # interaction per queue loop
            break

        # assuming next_inter still points to an interaction,
        # if any_found is False, readd to the queue (next_inter has been incr)
        # this is how we query all interactions one by one, even if we can't
        # fill it all immediately
        if not any_found and cmatch.next_inter < len(cmatch.interactions):
            queue.append(cmatch)
        elif not any_found:
            # if next_inter is already the last one, check if it is complete,
            # acceptable and non duplicate. If all three add to result list
            if cmatch.is_complete() and cmatch.is_acceptable():
                duplicate = False
                for res in results:
                    if cmatch.is_equal(res):
                        duplicate = True
                        break
                if not duplicate:
                    results.append(cmatch)

    return results
