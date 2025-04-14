from .mol_fragment import MolFragment


class ReactionTemplate:
    name: str
    reactants: list[str]
    distance_max: list[tuple[int, int, int, int, float]]
    distance_min: list[tuple[int, int, int, int, float]]
    angle_limits: list[tuple[int, int, int, int, int, int, float, float]]
    dihedral_limits: list[tuple[int, int, int, int, int, int, int, int, float, float]]
    probability: float
    global_counter: int
    global_limit: int
    break_groups: list[list[tuple[int, int]]]
    update_groups: list[list[tuple[int, int]]]
    product: MolFragment
    renames: list[tuple[int, int, str]]
    retypes: list[tuple[int, int, str]]
    recharges: list[tuple[int, int, float]]
    remasses: list[tuple[int, int, float]]

    def __init__(self, name):
        self.name = name
        self.reactants = []
        self.distance_max = []
        self.distance_min = []
        self.angle_limits = []
        self.dihedral_limits = []
        self.probability = 1.0
        self.global_counter = 0
        self.global_limit = None
        self.break_groups = []
        self.update_groups = []
        self.product = MolFragment(name)
        self.product.index_type = "pair"
        self.renames = []
        self.retypes = []
        self.recharges = []
        self.remasses = []

    def is_complete(self) -> bool:
        # when a reaction is finished parsing, if this returns False
        # it is considered an error

        return (
            len(self.reactants) > 0
        )
