from dataclasses import dataclass

@dataclass
class MoleculeType:
    molecule_name: str

    # type, res num, res name, atomname, charge, mass
    atoms: list[tuple[str, int, str, str, float | None, float | None]]
    exclusions: set[tuple[int, int]]
    interactions: list[tuple[str, list[int], list[float]]]

    # methods called by [molecules] and reactions
    def add_atoms_to_system(self, system) -> list[int]:
        res = []
        last_res = -1
        for type_, res_num, res_name, atom_name, charge, mass in self.atoms:
            if res_num != last_res:
                last_res = res_num
                system.new_residue()
            res.append(system.add_atom(
                atom_name, res_name, type_, charge, mass
            ))
        return res


    def instantiate(self, system, atom_indices: list[int]):
        # lookup force by name in system
        # call add interaction
        for (i, j) in self.exclusions.copy():
            if j > i:
                self.exclusions.add((j, i))
        for (i, j) in self.exclusions:
            if i > j:
                if i < 0 or j < 0:
                    # during reactions, missing optional atoms can do this
                    continue
                system.add_interaction(
                    "exclusion",
                    [atom_indices[i], atom_indices[j]],
                    []
                )
        for (name, members, params) in self.interactions:
            if any(x < 0 for x in members):
                # during reactions, missing optional atoms can do this
                continue
            system.add_interaction(
                name,
                [atom_indices[x] for x in members],
                params
            )