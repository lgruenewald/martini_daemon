from .reporter import Reporter
from ..utils import backup_try
import sys


def dump_topstar(topstar, file=sys.stdout):
    print("==== TopStar / Fragment Types ====", file=file)
    for k, molfrag in topstar.type_lookup.items():
        file.write(f"{k} ")
    file.write("\n")
    print("==== TopStar / GraphFragments ====", file=file)
    for g in topstar.graph_fragment_list:
        file.write(f"{g.name}: ({[name for (name, _, _, _) in g.atoms]}) ")
    file.write("\n")
    print("==== TopStar / ReactionTemplates ====", file=file)
    for rx in topstar.reaction_list:
        print(f"rx {rx.name} reactants {rx.reactants} products {rx.products}",
              file=file)
    print("==== TopStar / Fragments ====", file=file)
    for id, frag in topstar.frag_list.items():
        print(f"{id}: <frag {frag.name} ps {frag.particles} opt {frag.opt}>", file=file)
    print("==== TopStar / defrag list ====", file=file)
    for id, defrag in enumerate(topstar.defrag_list):
        print(f"particle {id} is in fragments {defrag}", file=file)
        if id > 100:
            break
#    print("==== TopStar / Interaction list ====", file=file)
#    for id, inter in enumerate(topstar.interaction_list):
#        print(f"particle {id} is in interactions {inter}", file=file)


class TopStarLogger(Reporter):
    """A reporter that dumps the state of T* after every modification algorithm
    run to topstar.log. Useful for debugging.
    """

    def init(self):
        backup_try("topstar.log")
        self.post_modification(-1)

    def post_modification(self, i):
        with open("topstar.log", "a") as file:
            print(f"===== Frame {i} =====", file=file)
            dump_topstar(self._topstar, file)


class ReactionReporter(Reporter):
    """A reporter that reports all reactions to reactions.log"""

    def init(self):
        backup_try("reactions.log")

    def pre_modification(self, reactions, i):
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

    def pre_detection(self, i):
        with open("fragment_counts.log", "a") as file:
            init_map = self._topstar.get_init_map()
            counts = [f"{key}: {len(values)}" for key, values in init_map.items()]
            print(f"Frame {i} {counts}", file=file)
