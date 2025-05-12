from .reporter import Reporter
from ..utils import backup_try
import sys
from ..topstar import TopStar


def dump_topstar(topstar: TopStar, file=sys.stdout) -> None:
    print("==== TopStar / Fragment Types ====", file=file)
    for k, molfrag in topstar.type_lookup.items():
        file.write(f"{k} ")
    file.write("\n")
    print("==== TopStar / GraphFragments ====", file=file)
    for g in topstar.graph_fragment_map.values():
        file.write(f"{g.name}: ({[name for (name, _, _, _) in g.atoms]}) ")
    file.write("\n")
    print("==== TopStar / ReactionTemplates ====", file=file)
    for _, rl in topstar.reactions.items():
        for rx in rl:
            print(f"rx {rx.name} reactants {rx.reactants}", file=file)
    print("==== TopStar / Fragments ====", file=file)
    for id, frag in topstar.frag_list.items():
        print(f"{id}: <frag {frag.name} ps {frag.atoms}>", file=file)
    print("==== TopStar / defrag list ====", file=file)
    for id, defrag in enumerate(topstar.defrag_list):
        print(f"particle {id} is in fragments {defrag}", file=file)
        if id > 100:
            break


class TopStarLogger(Reporter):
    """A reporter that dumps the state of T* after every modification algorithm
    run to <name>.toplog. Useful for debugging.
    """

    def init_dm(self, name):
        backup_try(f"{name}.toplog")

    def pre_detection(self, i, name) -> None:
        if i == 0:
            self.post_modification(i, name)

    def post_modification(self, i, name) -> None:
        with open(f"{name}.toplog", "a") as file:
            print(f"===== Frame {i} =====", file=file)
            dump_topstar(self._topstar, file)


class ReactionReporter(Reporter):
    """A reporter that reports all reactions to <name>.reactions"""

    def init_dm(self, name):
        backup_try(f"{name}.reactions")

    def pre_modification(self, i, reactions, name) -> None:
        with open(f"{name}.reactions", "a") as file:
            print(f"Frame {i}", file=file)
            for (frags, rx) in reactions:
                frags = [(frag.name, frag.frag_id, frag.atoms) for frag in frags]
                print(
                    f"Reaction {rx.name} reactants {frags}",
                    file=file
                )

    def interactive_line(self) -> str:
        return f"reactions: {self._simulation.reactions}"


class FragCountReporter(Reporter):
    """A reporter that logs the number of all fragments in T* at a given time
    to <name>.frags"""

    def init_dm(self, name):
        backup_try(f"{name}.frags")

    def pre_detection(self, i, name) -> None:
        data = ",".join(
            [f"{k}:{v}" for k, v in self._topstar.frag_counts.items()]
        )
        with open(f"{name}.frags", "a") as file:
            print(f"Frame:{i},{data}", file=file)

    def interactive_line(self) -> str:
        return f"fragments: {len(self._topstar.frag_list)}"
