from ..__core import MoleculeType, System
from ..__parser import ParseException


class ModificationTemplate(MoleculeType):
    def __init__(self) -> None:
        super().__init__()
        # reactants directive validates that the graphs exist btw
        self.reactants = []
        # flattened indices, so instantiate actually does things nicely
        self.break_groups: list[list[int]] = []
        self.update_groups: list[list[int]] = []
        self.renames: list[tuple[int, str]] = []
        self.retypes: list[tuple[int, str]] = []
        self.recharges: list[tuple[int, float]] = []
        self.remasses: list[tuple[int, float]] = []
        self.soft_core: list[tuple[int, float, float]] = []

    def add_atoms_to_system(self, system: System) -> None:
        """Add atoms to the system.

        Invalid to call for reactions, therefore it will raise an exception.
        """
        # reactions don't add atoms to system anymore
        raise ParseException(
            "Attempt to add a reaction to a [system]/[molecules]. Please use a [moleculetype] molecule."
        )

    def process_nrexcl(self) -> None:
        """Process nr_excl for a reaction.

        Only self.nrexcl == 1 is allowed for modification templates.
        """
        assert self.nrexcl == 1
        super().process_nrexcl()

    def instantiate(self, system: System, atom_indices: list[int]) -> None:
        """Apply a modification template on selected (flattened) atom indices."""
        # breaking bonds
        for break_group in self.break_groups:
            system.break_group(
                {atom_indices[x] for x in break_group if atom_indices[x] >= 0}
            )
        for update_group in self.update_groups:
            system.update_group(
                {atom_indices[x] for x in update_group if atom_indices[x] >= 0}
            )

        # [bonds], [angles]...
        super().instantiate(system, atom_indices)

        # changing atoms
        for i, name in self.renames:
            if atom_indices[i] >= 0:
                system.rename(atom_indices[i], name)
        for i, type_ in self.retypes:
            if atom_indices[i] >= 0:
                system.retype(atom_indices[i], type_)
        for i, charge in self.recharges:
            if atom_indices[i] >= 0:
                system.recharge(atom_indices[i], charge)
        for i, mass in self.remasses:
            if atom_indices[i] >= 0:
                system.remass(atom_indices[i], mass)

    def toggle_softcore(self, system: System, atom_indices: list[int], on: bool) -> None:
        """Toggle soft core on/off for minimization based on the reacting atom_indices."""
        for i, sc_lam, sc_alpha in self.soft_core:
            if atom_indices[i] >= 0:
                if on:
                    system.update_sc(atom_indices[i], sc_lam, sc_alpha)
                else:
                    system.update_sc(atom_indices[i], 1.0, 0.5)
