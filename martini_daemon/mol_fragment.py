from .forces.force import Force

class MolFragment:
    molecule_name: str
    # atoms: type, resnum, resname, atomname, chargegr, charge, mass
    atoms: list[tuple[str, int, str, str, int, float, float]]
    # exclusions
    exclusions: set[tuple[int | tuple[int, int], int | tuple[int, int]]]
    # interactions: generic members and params
    interactions: list[tuple[Force, list[int | tuple[int, int]], list[float]]]
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

    def save(self, f) -> None:
        f.dump(self.molecule_name)
        f.dump(self.atoms)
        f.dump(self.exclusions)
        f.dump(self.index_type)

        f.dump(len(self.interactions))
        for force, indices, params in self.interactions:
            f.dump(force.get_index())
            f.dump(force.get_class_name())
            f.dump(indices)
            f.dump(params)

    def load(self, f, sys) -> None:
        self.molecule_name = f.load()
        self.atoms = f.load()
        self.exclusions = f.load()
        self.index_type = f.load()

        for i in range(f.load()):
            force = sys.modular_forces[f.load()]
            assert f.load() == force.get_class_name(), "Attempt to load a checkpoint from a different version of daemon."
            indices = f.load()
            params = f.load()
            self.interactions.append((force, indices, params))
