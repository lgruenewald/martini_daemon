from .directive import Directive
from .gromacs_top_file import GromacsTopFile, register_directive
from .token_list import TokenList, TokenParseException


@register_directive
class NonbondParams(Directive):
    def line(self, tokens: TokenList) -> None:
        assert isinstance(self.parent, GromacsTopFile)
        type1 = tokens.unwrap(0, "word")
        type2 = tokens.unwrap(1, "word")
        funct = tokens.unwrap(2, "int")
        if funct not in {1}:
            raise TokenParseException(tokens[2], f"Unsupported function type {type}.")
        sigma = tokens.unwrap(3, "float")
        epsilon = tokens.unwrap(4, "float")
        if (nb_types := self.parent.system.additional_data.get("nb_types")) is None:
            nb_types = {}
            self.parent.system.additional_data["nb_types"] = nb_types
        nb_types[(type1, type2)] = (sigma, epsilon)

    def finish(self) -> None:
        pass

    @classmethod
    def is_mandatory(cls) -> bool:
        return True

    @classmethod
    def is_unique(cls) -> bool:
        return False

    @classmethod
    def is_valid_parent(cls, parent: Directive) -> bool:
        return type(parent) is GromacsTopFile

    @classmethod
    def get_name(cls) -> str:
        return "nonbond_params"
