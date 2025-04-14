from .graph import GraphFragment

class Fragment():
    name: str
    atoms: list[int]
    frag_id: int
    graph: GraphFragment

    def __init__(self, graph, id: int):
        self.name = graph.name
        self.atoms = []  # missing optional -> -1
        self.frag_id = id
        self.graph = graph

    def index_atom(self, i: int) -> int:
        return self.atoms[i]