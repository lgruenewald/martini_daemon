class ReactionTemplate:
    name: str
    reactants: list[str]
    distance_max: list[tuple[int, int, int, int, float]]
    distance_min: list[tuple[int, int, int, int, float]]
    angle_limits: list[tuple[int, int, int, int, int, int, float, float]]
    dihedral_limits: list[tuple[int, int, int, int, int, int, int, int, float, float]]
    relative_rate: float
    reaction_counter: int
    observed_rate: float | None
    break_groups: list[list[tuple[int, int]]]
    update_groups: list[list[tuple[int, int]]]
    renames: list[tuple[int, int, str]]
    retypes: list[tuple[int, int, str]]
    recharges: list[tuple[int, int, float]]
    remasses: list[tuple[int, int, float]]
    soft_core: list[tuple[int, int, float, float]]

    def __init__(self, name):
        self.name = name
        self.reactants = []
        self.distance_max = []
        self.distance_min = []
        self.angle_limits = []
        self.dihedral_limits = []
        self.reaction_counter = 0
        self.observed_rate = None
        self.relative_rate = None
        self.break_groups = []
        self.update_groups = []
        self.renames = []
        self.retypes = []
        self.recharges = []
        self.remasses = []
        self.soft_core = []

    def is_complete(self) -> bool:
        # when a reaction is finished parsing, if this returns False
        # it is considered an error

        return (
            len(self.reactants) > 0
        )
