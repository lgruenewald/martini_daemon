from ..__rust import TokenList, TokenParseException
from .directive import Directive
from .gromacs_top_file import register_directive
from .molecule_type_directive import MoleculeTypeDirective


@register_directive
class AtomsDirective(Directive):
    def line(self, tokens: TokenList) -> None:
        assert isinstance(self.parent, MoleculeTypeDirective)
        index = tokens.unwrap(0, "index")
        type_ = tokens.unwrap(1, "word")
        resid = tokens.unwrap(2, "word")
        # FIXME: better resid handling -> don't just drop the A
        if resid[-1] == "A":
            resid = resid[:-1]
        try:
            res_num = int(resid)
        except ValueError:
            raise TokenParseException(
                tokens[2],
                f"Residue ID should be a number, optionally followed by 'A'. Got: '{resid}'.",
            )
        res_name = tokens.unwrap(3, "word")
        atom_name = tokens.unwrap(4, "word")
        _ = tokens.unwrap(5, "int")  # charge group number
        charge = tokens.unwrap(6, "float", False, None)
        mass = tokens.unwrap(7, "float", False, None)
        if index != len(self.parent.molecule_type.atoms):
            raise TokenParseException(
                tokens[0],
                "Bad atom ID, are they out of order?"
                f" got id {index} but expected {len(self.parent.molecule_type.atoms)}",
            )
        self.parent.molecule_type.atoms.append(
            (type_, res_num, res_name, atom_name, charge, mass)
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
        return "atoms"
