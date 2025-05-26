from .reporter import Reporter


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
        self._print("==== TopStar / Fragment Types ====")
        for k, molfrag in topstar.type_lookup.items():
            self._write(f"{k} ")
        self._write("\n")
        self._print("==== TopStar / GraphFragments ====")
        for g in topstar.graph_fragment_map.values():
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
            self._print(f"particle {id} is in fragments {defrag}")
            if id > 100:
                break


class ReactionReporter(Reporter):
    """A reporter that reports all reactions to <name>.reactions"""

    def init_dm(self, name):
        self._open(name + ".reactions")

    def pre_modification(self, i, reactions, name) -> None:
        print(f"Frame {i}", file=self._handle)
        for (frags, rx) in reactions:
            frags = [(frag.name, frag.frag_id, frag.atoms) for frag in frags]
            self._print(
                f"Reaction {rx.name} reactants {frags}"
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
