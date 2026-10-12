from ..__rust import TokenList, TokenParseException
from .directive import Directive
from .gromacs_top_file import GromacsTopFile, register_directive


@register_directive
class DihedralTypesDirective(Directive):
    def line(self, tokens: TokenList) -> None:
        a = tokens.unwrap(0, "word")
        b = tokens.unwrap(1, "word")
        c = tokens.unwrap(2, "word")
        d = tokens.unwrap(3, "word")
        type_ = tokens.unwrap(4, "int")
        params = [tokens[i] for i in range(5, len(tokens))]
        assert isinstance(self.parent, GromacsTopFile)
        try:
            self.parent.system.add_dihedral_type(a, b, c, d, type_, params)
        except ValueError as e:
            raise TokenParseException(tokens[0], f"[{self.get_name()}]: {e}") from None

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
        return type(parent) is GromacsTopFile

    @classmethod
    def get_name(cls) -> str:
        return "dihedraltypes"
