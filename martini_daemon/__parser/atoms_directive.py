from typing import Any

from .directive import Directive
from .gromacs_top_file import register_directive
from .molecule_type_directive import MoleculeTypeDirective
from .token_list import TokenList, TokenParseException


@register_directive
class AtomsDirective(Directive):
    def __init__(self, parent: MoleculeTypeDirective, path: str, line_num: int) -> None:
        self.parent = parent
        self.path = path
        self.line_num = line_num

    def where(self) -> tuple[str, int]:
        return self.path, self.line_num

    def line(self, tokens: TokenList) -> None:
        index = tokens.unwrap(0, "index")
        type_ = tokens.unwrap(1, "word")
        res_num = tokens.unwrap(2, "int")
        res_name = tokens.unwrap(3, "word")
        atom_name = tokens.unwrap(4, "word")
        _ = tokens.unwrap(5, "int") # charge group number
        charge = tokens.unwrap(6, "float", None)
        mass = tokens.unwrap(7, "float", None)
        if index != len(self.parent.atoms):
            raise TokenParseException(
                tokens[0],
                "Bad atom ID, are they out of order?"
                f" got id {index} but expected {len(self.parent.atoms)}"
            )
        self.parent.atoms.append(
            (type_, res_num, res_name, atom_name, charge, mass)
        )


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
        return type(parent) is MoleculeTypeDirective

    @classmethod
    def get_name(cls) -> str:
        return "atoms"