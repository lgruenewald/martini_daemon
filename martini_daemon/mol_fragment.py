from .forces.force import Force, Interaction

class MolFragment:
    molecule_name: str
    # atoms: type, resnum, resname, atomname, chargegr, charge, mass
    atoms: list[tuple[str, int, str, str, int, float, float]]
    # tuple[int, str] instead of int when reaction TODO
    exclusions: set[tuple[int, int]]
    # interactions: generic members and params
    # tuple[int, str] instead of int when reaction
    interactions: list[tuple[Force, list[int], list[float]]]
    index_type = "index"

    def __init__(self, name):
        self.molecule_name = name
        self.atoms = []
        self.exclusions = set()
        self.interactions = []
        self.renames = []
        self.retypes = []

    def add_exclusion(self, i, j) -> None:
        """Adds an exclusion to the list of exclusions
        Does not add duplicate exclusions."""

        self.exclusions.add((i, j))
        self.exclusions.add((j, i))
