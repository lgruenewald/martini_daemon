from .old_graph import Graph


class Fragment():
    # TODO if detection2.py works out, change this to already
    # be a condensed version to reduce data massaging in D/M
    name: str
    atoms: list[int]
    frag_id: int
    graph: Graph

    def __init__(self, graph, id: int):
        self.name = graph.name
        self.atoms = []  # missing optional -> -1
        self.frag_id = id
        self.graph = graph

    def index_atom(self, i: int) -> int:
        return self.atoms[i]
