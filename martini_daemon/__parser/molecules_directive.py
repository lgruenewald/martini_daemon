from .directive import Directive
from .gromacs_top_file import GromacsTopFile, register_directive
from .token_list import TokenList


@register_directive
class MoleculesDirective(Directive):
    def line(self, tokens: TokenList) -> None:
        assert isinstance(self.parent, GromacsTopFile)
        self.parent.system.initial_molecules.append(
            (tokens.unwrap(0, "word"), tokens.unwrap(1, "int"))
        )

    def finish(self) -> None:
        pass

    @classmethod
    def is_mandatory(cls) -> bool:
        return True

    @classmethod
    def is_unique(cls) -> bool:
        return True

    @classmethod
    def is_valid_parent(cls, parent: Directive) -> bool:
        return type(parent) is GromacsTopFile

    @classmethod
    def get_name(cls) -> str:
        return "molecules"
