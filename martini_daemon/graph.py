# helper function for T* graph based fragments
from __future__ import annotations
from enum import Enum
from .forces.force import Interaction
from .sysstar import SysStar
from fnmatch import fnmatch


class GraphAtomType(Enum):
    NORMAL = 0
    OPT = 1
    NOT = 2


class GraphMatch():
    # mapping of atoms -> part_id
    graph: GraphFragment
    atoms: dict[str, int]
    # interaction matches, lists of None or Interaction
    interactions: list[None | Interaction]
    # reverse of atoms
    rev_atoms: dict[int, str]
    matched_inter: set[Interaction]

    def __init__(self, graph: GraphFragment):
        self.graph = graph
        self.atoms = {}
        self.interactions = [None for _ in graph.interactions]
        self.rev_atoms = {}
        self.matched_inter = set()

    def copy(self) -> GraphMatch:
        res = GraphMatch(self.graph)
        res.atoms = self.atoms.copy()
        res.interactions = self.interactions.copy()
        res.rev_atoms = self.rev_atoms.copy()
        res.matched_inter = self.matched_inter.copy()
        return res

    def add_atom(self, name: str, part_num: int):
        assert self.atoms.get(name) is None and part_num not in self.rev_atoms
        self.atoms[name] = part_num
        self.rev_atoms[part_num] = name

        # TODO -> have a neighbor list in GraphMatch and update it when add_atom
        # TODO -> convert interactions involved to a graph, keep it updated when
        # adding atoms

    def add_inter(self, id: int, inter: Interaction):
        assert self.interactions[id] is None and inter not in self.matched_inter
        self.interactions[id] = inter
        self.matched_inter.add(inter)

    def is_complete(self):
        # interactions are checked during construction
        # only partials that this should be called on are ones where
        # every time a new atom is added, all interactions are enforced
        # that can be enforced when adding that atom // TLDR only checks atoms
        for (name, _, _, part_type) in self.graph.atoms:
            found = self.atoms.get(name) is not None
            if part_type == GraphAtomType.NORMAL and not found:
                return False
        return True

    def is_acceptable(self):
        # is complete + not checking
        for (name, _, _, part_type) in self.graph.atoms:
            found = self.atoms.get(name) is not None
            if part_type == GraphAtomType.NORMAL and not found:
                return False
            if part_type == GraphAtomType.NOT and found:
                return False
        return True

    def is_equal(self, other: GraphMatch) -> bool:
        # also see note for is_acceptable -> only checks atoms
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


class GraphFragment():
    name: str
    # only summon this graph on atoms of this molecule during its instantion
    # and only during the start of the simulation (TODO make this more flexible)
    # if empty, all molecules and during reactions too
    molecules: list[str]
    # list[(part_id, name_pat, type_pat, type)]
    atoms: list[tuple[str, str, str, GraphAtomType]]
    # list[(interaction_type, list[part_id])]
    interactions: list[tuple[str, list[str]]]
    equivalents: list[set[str]]

    def __init__(self, name):
        self.name = name
        self.molecules = []
        self.atoms = []
        self.interactions = []
        self.equivalents = []

    def get_neighbors(self, part: int, interactions: list[list[Interaction]]):
        neighbors = set()
        for inter in interactions[part]:
            for member in inter.get_members():
                neighbors.add(member)
        return neighbors

    def check_atom_interactions(self, part_name: str, partial: GraphMatch,
                                interactions: list[list[Interaction]]):
        # TODO accelerate this using the partial's built in interaction graph
        """
            Goes through all interactions of part_name in the partial match.
            Adds all interactions that it can to partial (mutating it).

            Returns whether all interactions that should be fulfilled by
            now (between particles already in partial) have been fulfilled.
            (when False this partial can be considered a bad match, so it
            returns early)
        """
        part_num = partial.atoms[part_name]
        inters = interactions[part_num]
        # enumerate over interactions in GraphFragment
        for i, (inter_type, inter_parts) in enumerate(self.interactions):
            # only consider unfulfilled interactions in graph (PartialMatch)
            if partial.interactions[i] is not None:
                continue
            # only consider interactions containing the particle
            if part_name not in inter_parts:
                continue
            # if there is a missing member, then it's too early to check
            part_nums = [partial.atoms.get(ipart) for ipart in inter_parts]
            if any(map(lambda x: x is None, part_nums)):
                continue
            part_nums = set(part_nums)
            found = False
            for inter in inters:
                if inter in partial.matched_inter:
                    continue
                # type filter
                if not inter.is_instance(inter_type):
                    continue
                # members filter
                inter_part_nums = set(inter.get_members())
                if part_nums != inter_part_nums:
                    continue
                # got here? it's a match so add it
                found = True
                partial.add_inter(i, inter)
                break
            if not found:
                return False
        return True

    def try_match(self, particles: set[int], partial: GraphMatch,
                  sysstar: SysStar, interactions: list[list[Interaction]],
                  ) -> list[GraphMatch]:
        # TODO rewrite documentation strings
        # TODO interaction graph informed new neighbor picking
        """
            Please use match_particles() from outside

            Warning: it will also return not acceptable matches containing NOT
            atoms. Filter it with .is_acceptable().

            Return all good fragment matches that can be built up from
            particles and partial, that do not exist in the fragment
            list yet.

            Works by trying to add a single particle to partial at a time
            and recursively calling itself.

            This is a greedy algorithm, that is it will match as many atoms
            as it can and not return partial matches. This is useful for
            properly supporting optionals and nots.

            The atoms in the graph have no defined order. Matching optionals
            before nots is not a necessity, since in case of ambiguity it should
            produce multiple graphs with all combinations, and one of them
            would have the atoms arranged in a way to prioritize optionals.

        """

        # list of stuff to recurse on later
        partials = []
        new_particles = []
        part_nums = []
        # particles to still check in all children
        keep = set()

        # filter through particles and collect new recursions
        for part_num in particles:
            if part_num in partial.rev_atoms.keys():
                continue
            name, type = sysstar.get_particle_name_type(part_num)
            neighbors = self.get_neighbors(part_num, interactions)
            for (part_name, name_pat, type_pat, _) in self.atoms:
                # ignore graph atoms already matched
                if partial.atoms.get(part_name) is not None:
                    continue
                # don't match wrong name/type
                if not fnmatch(name, name_pat) or not fnmatch(type, type_pat):
                    continue
                # unfulfilled atom with a promising name and type
                cpartial = partial.copy()
                cpartial.add_atom(part_name, part_num)
                if not self.check_atom_interactions(part_name, cpartial,
                                                    interactions):
                    continue

                # this algorithm always finishes because during each
                # recursion we must always add one atom to the matched list

                # this particle can stay in particles if we can still grow the
                # graph this way in one of the recursions
                keep.add(part_num)
                partials.append(cpartial)
                # beware, multiple copies of the same set
                new_particles.append(neighbors)
                part_nums.append(part_num)

        res = []
        for i in range(len(partials)):
            cpartial = partials[i]
            cparticles = (keep | new_particles[i]) - {part_nums[i]}
            new_matches = self.try_match(cparticles, cpartial, sysstar,
                                         interactions)
            # don't add duplicates
            for m in new_matches:
                if all(map(lambda x: not m.is_equal(x), res)):
                    res.append(m)

        # greedy algorithm => only check if the current graph match is
        # acceptable when there isn't a match accepted that this is a subset of
        if len(res) == 0:
            if partial.is_complete():
                res.append(partial)

        return res

    def match_particles(self, particles: set[int],
                        sysstar: SysStar,
                        interactions: list[list[Interaction]],
                        ) -> list[GraphMatch]:
        # TODO find starting matches of a single atom separately here
        # then call try_match on all
        # TODO remove duplicates here
        """
            Friendly wrapper around try_match + filter out results with NOT

            see try_match
        """
        partial = GraphMatch(self)
        matches = self.try_match(particles, partial, sysstar, interactions)
        return list(filter(lambda x: x.is_acceptable(), matches))
