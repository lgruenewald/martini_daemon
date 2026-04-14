from typing import Any


class MoleculeType:
    def __init__(self):
        """A class that contains all the information that is in a moleculetype.
        Along with helpers to mutate system to add the required atoms for it
        to system, and then to instantiate the bonded interactions on top
        of the new atoms added.
        """
        self.name: str | None = None
        self.nrexcl: int | None = None

        # type, res num, res name, atomname, charge, mass
        self.atoms: list[tuple[str, int, str, str, float | None, float | None]] = []
        self.exclusions: set[tuple[int, int]] = set()
        # name, members, params, generates_excl?
        self.interactions: list[tuple[str, list[int], list[float], bool]] = []

    # methods called by [molecules] and reactions
    def add_atoms_to_system(self, system: Any) -> list[int]:
        """Adds the atoms in molecule type to system.

        :param system: The Martini Daemon system.
        :return: list of atom indices that got added to the system.
        """
        res = []
        last_res = -1
        for type_, res_num, res_name, atom_name, charge, mass in self.atoms:
            if res_num != last_res:
                last_res = res_num
                system.new_residue()
            res.append(system.add_atom(atom_name, res_name, type_, charge, mass))
        return res

    def process_nrexcl(self) -> None:
        """Process bond-generated exclusions.

        Called after parsing. It looks at self.interactions to see which bonds generate exclusions, and then
        it adds the generated exclusions to it based on the nr_excl argument in [moleculetype].
        """
        for _, members, _, is_excl in self.interactions:
            if is_excl:
                assert len(members) == 2
                self.exclusions.add((members[0], members[1]))

    def instantiate(self, system, atom_indices: list[int]) -> None:
        """Adds the bonded interactions and exclusions stored in this MoleculeType to system.

        :param system: The Martini Daemon system.
        :param atom_indices: the selected atom indices belonging to this molecule.
        """
        # lookup force by name in system
        # call add interaction
        for i, j in self.exclusions.copy():
            if j > i:
                self.exclusions.add((j, i))
        for i, j in self.exclusions:
            if i > j:
                if i < 0 or j < 0:
                    # during reactions, missing optional atoms can do this
                    continue
                system.add_interaction(
                    "exclusion", [atom_indices[i], atom_indices[j]], []
                )
        for name, members, params, _ in self.interactions:
            if any(atom_indices[x] < 0 for x in members):
                # during reactions, missing optional atoms can do this
                continue
            system.add_interaction(name, [atom_indices[x] for x in members], params)
