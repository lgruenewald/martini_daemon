from typing import Any

from .parser import ParseException
from .directive import  Directive
from .token_list import TokenList, TokenParseException
from .gromacs_top_file import GromacsTopFile, register_directive
from ..__core import MoleculeType

@register_directive
class MoleculeTypeDirective(Directive):
    def __init__(self, parent: GromacsTopFile, path: str, line_num: int) -> None:
        super().__init__(parent, path, line_num)

        self.molecule_name: str | None = None
        self.nrexcl: int | None = None
        self.system = parent.system

        # type, res num, res name, atomname, charge, mass
        self.atoms: list[tuple[str, int, str, str, float | None, float | None]] = []
        self.exclusions: set[tuple[int, int]] = set()
        self.interactions: list[tuple[str, list[int], list[float], bool]] = []

    def line(self, tokens: TokenList) -> None:
        self.molecule_name = tokens.unwrap(0, "word")
        self.nrexcl = tokens.unwrap(1, "int")
        if self.nrexcl < 1 or self.nrexcl > 4:
            raise TokenParseException(tokens[1], "Only nr. excl. between 1 and 4 allowed.")
        tokens.assert_no_more_than(2)

    def finish(self):
        if self.molecule_name is None or self.nrexcl is None:
            raise ParseException("Molecule name and nr. excl. expected.")
        self.process_nrexcl()
        self.parent.system.molecule_types[self.molecule_name] = MoleculeType(
            self.molecule_name,
            self.atoms,
            self.exclusions,
            [(force, members, params) for force, members, params, _ in self.interactions]
        )

    @classmethod
    def is_mandatory(cls):
        return False

    @classmethod
    def is_unique(cls):
        return False

    @classmethod
    def is_valid_parent(cls, parent: Any) -> bool:
        return isinstance(parent, GromacsTopFile)

    @classmethod
    def get_name(cls) -> str:
        return "moleculetype"

    def process_nrexcl(self):
        assert self.nrexcl == 1, "Martini uses nrexcl == 1" # TODO temporary
        # TODO make a nrexcl = 2 test / some AA force field tests for correctness
        for _, members, _, is_excl in self.interactions:
            if is_excl:
                assert len(members) == 2
                self.exclusions.add((members[0], members[1]))

    # methods called by [bonds], [angles]...
    def parse_index(self, tokens: TokenList, index: int) -> int:
        return tokens.unwrap(index, "index")