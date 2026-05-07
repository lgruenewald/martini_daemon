from typing import Any

from ..__core import MoleculeType
from .directive import Directive
from .gromacs_top_file import GromacsTopFile, register_directive
from .parser import ParseException
from .token_list import TokenList, TokenParseException


@register_directive
class MoleculeTypeDirective(Directive):
    def __init__(self, parent: GromacsTopFile, path: str, line_num: int) -> None:
        super().__init__(parent, path, line_num)
        self.molecule_type = MoleculeType()
        self.system = parent.system

    def line(self, tokens: TokenList) -> None:
        self.molecule_type.name = tokens.unwrap(0, "word")
        self.molecule_type.nrexcl = tokens.unwrap(1, "int")
        if self.molecule_type.nrexcl != 1:
            raise TokenParseException(tokens[1], "Only nr. excl. == 1 is supported.")
        tokens.assert_no_more_than(2)

    def finish(self) -> None:
        if self.molecule_type.name is None or self.molecule_type.nrexcl is None:
            raise ParseException("Molecule name and nr. excl. expected.")
        self.molecule_type.process_nrexcl()
        self.system.molecule_types[self.molecule_type.name] = self.molecule_type

    @classmethod
    def is_mandatory(cls) -> bool:
        return False

    @classmethod
    def is_unique(cls) -> bool:
        return False

    @classmethod
    def is_valid_parent(cls, parent: Any) -> bool:
        return isinstance(parent, GromacsTopFile)

    @classmethod
    def get_name(cls) -> str:
        return "moleculetype"

    # methods called by [bonds], [angles]...
    def parse_index(self, tokens: TokenList, index: int) -> int:
        return tokens.unwrap(index, "index")
