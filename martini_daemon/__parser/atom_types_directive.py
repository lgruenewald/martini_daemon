from .directive import Directive
from .gromacs_top_file import GromacsTopFile, register_directive
from .parser import ParseException
from .token_list import TokenList, TokenParseException


@register_directive
class AtomTypesDirective(Directive):
    def line(self, tokens: TokenList) -> None:
        if len(tokens) != 6:
            raise ParseException(
                "Only atomtypes lines formatted as type, m, q,"
                "atom_type, sigma, epsilon are supported."
            )
        type = tokens.unwrap(0, "word")
        mass = tokens.unwrap(1, "float")
        charge = tokens.unwrap(2, "float")
        atom_type = tokens.unwrap(3, "word")
        if atom_type != "A":
            raise TokenParseException(tokens[3], "Only 'A' atom type supported.")
        sigma = tokens.unwrap(4, "float")
        epsilon = tokens.unwrap(5, "float")
        if sigma != 0.0:
            raise TokenParseException(
                tokens[4], "Only sigma=0 is supported in [atomtypes]."
            )
        if epsilon != 0.0:
            raise TokenParseException(
                tokens[5], "Only epsilon=0 is supported in [atomtypes]."
            )
        self.parent.system.add_atom_type(type, charge, mass)

    def finish(self) -> None:
        pass

    @classmethod
    def is_mandatory(cls) -> bool:
        return True

    @classmethod
    def is_unique(cls) -> bool:
        return False

    @classmethod
    def is_valid_parent(cls, parent) -> bool:
        return type(parent) is GromacsTopFile

    @classmethod
    def get_name(cls) -> str:
        return "atomtypes"
