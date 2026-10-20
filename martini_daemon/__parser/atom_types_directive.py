from ..__rust import TokenList, TokenParseException
from .directive import Directive
from .gromacs_top_file import GromacsTopFile, register_directive
from .parser import ParseException


@register_directive
class AtomTypesDirective(Directive):
    def line(self, tokens: TokenList) -> None:
        assert isinstance(self.parent, GromacsTopFile)
        if len(tokens) == 6:
            type = tokens.unwrap(0, "word")
            mass = tokens.unwrap(1, "float")
            charge = tokens.unwrap(2, "float")
            atom_type = tokens.unwrap(3, "word")
            if atom_type != "A":
                raise TokenParseException(tokens[3], "Only 'A' atom type supported.")
            sigma = tokens.unwrap(4, "float")
            epsilon = tokens.unwrap(5, "float")
        elif len(tokens) == 7:
            type = tokens.unwrap(0, "word")
            _Z = tokens.unwrap(1, "float")
            mass = tokens.unwrap(2, "float")
            charge = tokens.unwrap(3, "float")
            atom_type = tokens.unwrap(4, "word")
            if atom_type != "A":
                raise TokenParseException(tokens[4], "Only 'A' atom type supported.")
            sigma = tokens.unwrap(5, "float")
            epsilon = tokens.unwrap(6, "float")
        # FIXME: bond_type column support
        else:
            raise ParseException(
                "Only atomtypes lines formatted as type, [Z], m, q,"
                "atom_type, sigma, epsilon are supported."
            )
        self.parent.system.add_atom_type(type, charge, mass, sigma, epsilon)

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
        return "atomtypes"
