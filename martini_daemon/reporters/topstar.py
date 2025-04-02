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
    for g in topstar.graph_fragment_list:
        file.write(f"{g.name}: ({[name for (name, _, _, _) in g.atoms]}) ")
    file.write("\n")
    print("==== TopStar / ReactionTemplates ====", file=file)
    for _, rl in topstar.reactions.items():
        for rx in rl:
            print(f"rx {rx.name} reactants {rx.reactants}", file=file)
    print("==== TopStar / Fragments ====", file=file)
    for id, frag in topstar.frag_list.items():
        print(f"{id}: <frag {frag.name} ps {frag.particles} opt {frag.opt}>", file=file)
    print("==== TopStar / defrag list ====", file=file)
    for id, defrag in enumerate(topstar.defrag_list):
        print(f"particle {id} is in fragments {defrag}", file=file)
        if id > 100:
            break


class TopStarLogger(Reporter):
    """A reporter that dumps the state of T* after every modification algorithm
    run to topstar.log. Useful for debugging.
    """

    def init(self):
        backup_try("topstar.log")
        self.post_modification(-1)

    def post_modification(self, i) -> None:
        with open("topstar.log", "a") as file:
            print(f"===== Frame {i} =====", file=file)
            dump_topstar(self._topstar, file)


class ReactionReporter(Reporter):
    """A reporter that reports all reactions to reactions.log"""

    def init(self):
        backup_try("reactions.log")

    def pre_modification(self, reactions, i) -> None:
        with open("reactions.log", "a") as file:
            print(f"Frame {i}", file=file)
            for (frags, rx) in reactions:
                frags = [(frag.name, frag.frag_id, frag.particles) for frag in frags]
                print(
                    f"Reaction {rx.name} reactants {frags}",
                    file=file
                )


class FragCountReporter(Reporter):
    """A reporter that logs the number of all fragments in T* at a given time
    to fragment_counts.log"""

    def init(self):
        backup_try("fragment_counts.log")

    def pre_detection(self, i) -> None:
        with open("fragment_counts.log", "a") as file:
            nums: dict[str, int] = {}
            for _, v in self._topstar.frag_list.items():
                name = v.name
                if nums.get(name) is None:
                    nums[name] = 1
                else:
                    nums[name] += 1

            counts = [f"{key}: {count}" for key, count in nums.items()]
            print(f"Frame {i} {counts}", file=file)
