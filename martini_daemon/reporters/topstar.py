from .reporter import Reporter
from ..helpers.monomer import generate_mapping


class TopStarLogger(Reporter):
    """A reporter that dumps the state of T* after every modification algorithm
    run to <name>.toplog. Useful for debugging.
    """

    def init_dm(self, name):
        self._open(name + ".toplog")

    def pre_detection(self, i, name) -> None:
        # first ever frame
        if i == 0:
            self.post_modification(i, name)

    def post_modification(self, i, name) -> None:
        topstar = self._topstar
        self._print(f"===== Frame {i} =====")
        self._print("==== TopStar / Molecules ====")
        for k, molfrag in topstar.molecules.items():
            self._write(f"{k} ")
        self._write("\n")
        self._print("==== TopStar / Graphs ====")
        for g in topstar.graphs.values():
            self._write(f"{g.name}: ({[name for (name, _, _, _) in g.atoms]}) ")
        self._write("\n")
        self._print("==== TopStar / ReactionTemplates ====")
        for _, rl in topstar.reactions.items():
            for rx in rl:
                self._print(f"rx {rx.name} reactants {rx.reactants}")
        self._print("==== TopStar / Fragments ====")
        for id, frag in topstar.frag_list.items():
            self._print(f"{id}: <frag {frag.name} ps {frag.atoms}>")
        self._print("==== TopStar / defrag list ====")
        for id, defrag in enumerate(topstar.defrag_list):
            self._print(f"atom {id} is in fragments {defrag}")
            if id > 100:
                break


class ReactionReporter(Reporter):
    """A reporter that reports all reactions to <name>.reactions"""

    def __init__(self, molid=False):
        """
        A reporter that reports all reactions to <name>.reactions.

        The created file has a text format, where every line is a reaction.
        First, the frame number and reaction name are separated by a comma,
        then, each reactant is separated by a semicolon. Each reactant
        will have its frag name, internal frag id, and atom indices
        (-1 for missing optional or forbidden atoms) printed.

        Example:
        frame,reaction_name;reactant1_name,reactant1_id,atoms...;...reactantn_name,reactantn_id,atoms...

        If molid is True, additionally the indices of initial molecules
        (see helpers/monomer) are printed in parentheses, prefixed with mol:
        after reactant IDs, before atoms.

        Example:
        frame,reaction_name;reactant1_name,reactant1_id(mol:molid1,...molidn),atoms...;...reactantn_name,reactantn_id(mol:molid1,...molidn),atoms...
        """
        self.molid = molid

    def init_dm(self, name):
        if self.molid:
            _, _, self.mapping = generate_mapping(self._topstar.initial_molecules)
        self._open(name + ".reactions")
        if not self.molid:
            self._print(
                "# frame,reaction_name;"
                "reactant1_name,reactant1_id,atoms...;..."
                "reactantn_name,reactantn_id,atoms..."
            )
        else:
            self._print(

                "# frame,reaction_name;"
                "reactant1_name,reactant1_id(mol:molid1,...molidn),atoms...;..."
                "reactantn_name,reactantn_id(mol:molid1,...molidn),atoms...;"
            )

    def get_molids(self, atoms):
        if not self.molid:
            return ""
        mols = set()
        for atom in atoms:
            if atom != -1:
                mols.add(f"{self.mapping[atom]}")
        return "(mol:" + ",".join(mols) + ")"

    def pre_modification(self, i, reactions, name) -> None:
        for (frags, rx) in reactions:
            self._print(
                f"{i},{rx.name};"
                + ";".join([
                    f"{frag.name},{frag.frag_id}"
                    f"{self.get_molids(frag.atoms)},"
                    + ",".join([
                        f"{atom}"
                        for atom in frag.atoms
                    ])
                    for frag in frags
                ])
            )

    def interactive_line(self) -> str:
        return f"reactions: {self._simulation.reactions}"


class FragCountReporter(Reporter):
    """A reporter that logs the number of all fragments in T* at a given time
    to <name>.frags"""

    def init_dm(self, name):
        self._open(name + ".frags")

    def pre_detection(self, i, name) -> None:
        data = ",".join(
            [f"{k}:{v}" for k, v in self._topstar.frag_counts.items()]
        )
        self._print(f"Frame:{i},{data}")

    def interactive_line(self) -> str:
        return f"fragments: {len(self._topstar.frag_list)}"


# for debugging reactions
class LastReactionGeometryReporter(Reporter):
    def pre_modification(self, i, reactions, name):
        self._sysstar.write_gro(name + "_last_reaction.gro")


class ReactionException(Exception):
    pass


# for debugging reactions, to stop simulations after reaction
class ReactionExceptionReporter(Reporter):
    def post_modification(self, i, name):
        raise ReactionException
