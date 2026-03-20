class ReactionTemplate:
    def __init__(self, name):
        self.name: str = name
        self.reactants: list[str] = []
        self.distance_max: list[tuple[int, int, int, int, float]] = []
        self.distance_min: list[tuple[int, int, int, int, float]] = []
        self.angle_limits: list[tuple[int, int, int, int, int, int, float, float]] = []
        self.dihedral_limits: list[tuple[int, int, int, int, int, int, int, int, float, float]] = []
        self.reaction_counter: int = 0
        self.observed_rate: None | float = None
        self.relative_rate: None | float = None
        self.break_groups: list[list[tuple[int, int]]] = []
        self.update_groups: list[list[tuple[int, int]]] = []
        self.renames: list[tuple[int, int, str]] = []
        self.retypes: list[tuple[int, int, str]] = []
        self.recharges: list[tuple[int, int, float]] = []
        self.remasses: list[tuple[int, int, float]] = []
        self.soft_core: list[tuple[int, int, float, float]] = []
        self.probability: None | float = None

    def is_complete(self) -> bool:
        # when a reaction is finished parsing, if this returns False
        # it is considered an error

        return (
            len(self.reactants) > 0
        )
