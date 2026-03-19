from typing import Any

from .directive import Directive
from .molecule_type_directive import MoleculeTypeDirective
from .gromacs_top_file import register_directive
from ..__parser import TokenList


@register_directive
class ExclusionsDirective(Directive):

    def line(self, tokens: TokenList) -> None:
        i = self.parent.parse_index(tokens, 0)
        # j is mandatory -> separate
        j = self.parent.parse_index(tokens, 1)
        self.parent.exclusions.add((i, j))
        for k in range(2, len(tokens)):
            self.parent.exclusions.add((
                i, self.parent.parse_index(tokens, k)
            ))


    def finish(self):
        pass

    @classmethod
    def is_mandatory(cls):
        return False

    @classmethod
    def is_unique(cls):
        return False

    @classmethod
    def is_valid_parent(cls, parent: Any) -> bool:
        return isinstance(parent, MoleculeTypeDirective)

    @classmethod
    def get_name(cls) -> str:
        return "exclusions"
