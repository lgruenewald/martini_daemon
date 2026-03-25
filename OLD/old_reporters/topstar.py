from .reporter import Reporter
from ..old_helpers.monomer import generate_mapping


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

    def post_modification(self, i, name, reactions) -> None:
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
                self._print(f"rx {rx.name} reactants {rx.__reactants}")
        self._print("==== TopStar / Fragments ====")
        for id, frag in topstar.frag_list.items():
            self._print(f"{id}: <frag {frag.name} ps {frag.atoms}>")
        self._print("==== TopStar / defrag list ====")
        for id, defrag in enumerate(topstar.defrag_list):
            self._print(f"atom {id} is in fragments {defrag}")
            if id > 100:
                break



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
