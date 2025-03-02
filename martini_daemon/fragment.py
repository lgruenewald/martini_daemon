class Fragment():
    name: str
    particles: list[int]
    frag_id: int

    def __init__(self, name, id):
        self.name = name
        self.particles = []
        self.frag_id = id
        self.graph = None
