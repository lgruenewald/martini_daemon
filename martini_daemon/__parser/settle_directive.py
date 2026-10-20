
from ..__rust import TokenList, TokenParseException
from .directive import Directive
from .gromacs_top_file import register_directive
from .molecule_type_directive import MoleculeTypeDirective


@register_directive
class SettlesDirective(Directive):
    def line(self, tokens: TokenList) -> None:
        atom = tokens.unwrap(0, "index")
        type_ = tokens.unwrap(1, "int")
        if type_ != 1:
            raise TokenParseException(
                tokens[1], f"Unsupported settles type: {type_}. Support types are: 1."
            )
        dOH = tokens.unwrap(2, "float")
        dHH = tokens.unwrap(3, "float")
        assert isinstance(self.parent, MoleculeTypeDirective)
        self.parent.molecule_type.interactions.append(
            ("constraint", [atom, atom + 1], [dOH], True)
        )
        self.parent.molecule_type.interactions.append(
            ("constraint", [atom, atom + 2], [dOH], True)
        )
        self.parent.molecule_type.interactions.append(
            ("constraint", [atom + 1, atom + 2], [dHH], True)
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
        return "settles"
