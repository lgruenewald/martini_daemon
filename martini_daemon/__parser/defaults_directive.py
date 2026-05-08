from .directive import Directive
from .gromacs_top_file import GromacsTopFile, register_directive
from .token_list import TokenList, TokenParseException


@register_directive
class Defaults(Directive):
    def line(self, tokens: TokenList) -> None:
        nb_type = tokens.unwrap(0, "int")
        if nb_type != 1:
            raise TokenParseException(
                tokens[0], f"Unsupported nonbonded type {nb_type}."
            )
        combination_rule = tokens.unwrap(1, "int")
        if combination_rule != 2:
            raise TokenParseException(
                tokens[1], f"Unsupported combination rule {combination_rule}."
            )
        tokens.assert_no_more_than(2)

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
        return "defaults"
