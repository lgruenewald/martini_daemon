from typing import Any

from .parser import ParseException
from .directive import  Directive
from .token_list import TokenList, TokenParseException
from .gromacs_top_file import GromacsTopFile


class MoleculeType(Directive):
    def __init__(self, parent: GromacsTopFile, path: str, line_num: int) -> None:
        self.path = path
        self.line_num = line_num
        self.molecule_name: str | None = None
        self.nrexcl: int | None = None
        self.parent = parent
        self.system = parent.system

        # type, res num, res name, atomname, charge_gr, charge, mass
        self.atoms: list[tuple[str, int, str, str, int, float | None, float | None]] = []
        self.exclusions: set[tuple[int, int]] = set()
        self.interactions: list[tuple[str, list[int], list[float], bool]] = []

    def where(self) -> tuple[str, int]:
        return self.path, self.line_num

    def line(self, tokens: TokenList) -> None:
        self.molecule_name = tokens.unwrap(0, "word")
        self.nrexcl = tokens.unwrap(1, "int")
        if self.nrexcl < 1 or self.nrexcl > 4:
            raise TokenParseException(tokens[1], "Only nr. excl. between 1 and 4 allowed.")
        tokens.assert_no_more_than(2)

    def finish(self):
        if self.molecule_name is None or self.nrel_excl is None:
            raise ParseException("Molecule name and nr. excl. expected.")
        self.process_nrexcl()
        self.parent.get("molecule_types", {})[self.molecule_name] = self

    @classmethod
    def is_mandatory(cls):
        return False

    @classmethod
    def is_unique(cls):
        return False

    @classmethod
    def is_valid_parent(cls, parent: Any) -> bool:
        return type(parent) is GromacsTopFile

    @classmethod
    def get_name(cls) -> str:
        return "moleculetype"

    def process_nrel_excl(self):
        raise NotImplementedError
    # TODO read out nrel excl from parent, TODO check spelling
    # go through interactions, call the Force API to whether exclusions should be generated
    # do the nrel excl algo, look up openmm
    # add exclusions to self
    # TODO make a nrel excl = 2 test / some AA force field tests for correctness

    # methods called by [bonds], [angles]...
    def parse_index(self, tokens: TokenList, index: int) -> int:
        return tokens.unwrap(index, "index")

    # methods called by [molecules] and reactions
    def add_atoms_to_system(self, system):
        raise NotImplementedError

    def instantiate(self, system, atom_indices: list[int]):
        # lookup force by name in system
        # call add interaction
        raise NotImplementedError