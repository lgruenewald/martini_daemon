from typing import Any

from .directive import Directive
from .gromacs_top_file import register_directive
from .token_list import TokenList
from .gromacs_top_file import GromacsTopFile


@register_directive
class MoleculesDirective(Directive):

    def line(self, tokens: TokenList) -> None:
        self.parent.system.initial_molecules.append((
            tokens.unwrap(0, "word"),
            tokens.unwrap(1, "int")
        ))

    def finish(self):
        pass

    @classmethod
    def is_mandatory(cls):
        return True

    @classmethod
    def is_unique(cls):
        return True

    @classmethod
    def is_valid_parent(cls, parent: Any) -> bool:
        return type(parent) == GromacsTopFile

    @classmethod
    def get_name(cls) -> str:
        return "molecules"