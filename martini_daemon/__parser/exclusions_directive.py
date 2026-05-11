from .directive import Directive
from .gromacs_top_file import register_directive
from .molecule_type_directive import MoleculeTypeDirective
from .token_list import TokenList


@register_directive
class ExclusionsDirective(Directive):
    def line(self, tokens: TokenList) -> None:
        assert isinstance(self.parent, MoleculeTypeDirective)
        i = self.parent.parse_index(tokens, 0)
        # j is mandatory -> separate
        j = self.parent.parse_index(tokens, 1)
        self.parent.molecule_type.exclusions.add((i, j))
        for k in range(2, len(tokens)):
            self.parent.molecule_type.exclusions.add(
                (i, self.parent.parse_index(tokens, k))
            )

    def finish(self) -> None:
        pass

    @classmethod
    def is_mandatory(cls) -> bool:
        return False

    @classmethod
    def is_unique(cls) -> bool:
        return False

    @classmethod
    def is_valid_parent(cls, parent: Directive) -> bool:
        return isinstance(parent, MoleculeTypeDirective)

    @classmethod
    def get_name(cls) -> str:
        return "exclusions"
