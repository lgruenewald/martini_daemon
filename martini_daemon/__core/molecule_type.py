from dataclasses import dataclass

@dataclass
class MoleculeType:
    molecule_name: str

    # type, res num, res name, atomname, charge_gr, charge, mass
    atoms: list[tuple[str, int, str, str, int, float | None, float | None]]
    exclusions: set[tuple[int, int]]
    interactions: list[tuple[str, list[int], list[float]]]

    # methods called by [molecules] and reactions
    def add_atoms_to_system(self, system):
        raise NotImplementedError

    def instantiate(self, system, atom_indices: list[int]):
        # lookup force by name in system
        # call add interaction
        raise NotImplementedError