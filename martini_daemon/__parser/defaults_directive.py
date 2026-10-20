from ..__rust import TokenList, TokenParseException
from .directive import Directive
from .gromacs_top_file import GromacsTopFile, register_directive


@register_directive
class Defaults(Directive):
    def line(self, tokens: TokenList) -> None:
        assert isinstance(self.parent, GromacsTopFile)
        if self.parent.system.additional_data.get("nb_type") is not None:
            raise TokenParseException(
                tokens[0], "[defaults] were already specified earlier."
            )
        nb_type = tokens.unwrap(0, "int")
        if nb_type != 1:
            raise TokenParseException(
                tokens[0], f"Unsupported nonbonded type {nb_type}."
            )
        combination_rule = tokens.unwrap(1, "int")
        if combination_rule not in {1, 2, 3}:
            raise TokenParseException(
                tokens[1],
                f"Unsupported combination rule {combination_rule}. Only comb. rules 1, 2 and 3 are supported.",
            )
        gen_pairs = tokens.unwrap(2, "word", mandatory=False, default="no").lower()
        if gen_pairs not in {"yes", "no"}:
            raise TokenParseException(tokens[2], "gen-pairs must be 'yes' or 'no'.")
        fudgeLJ = tokens.unwrap(3, "float", mandatory=False, default=1.0)
        fudgeQQ = tokens.unwrap(4, "float", mandatory=False, default=1.0)
        self.parent.system.additional_data["nb_type"] = nb_type
        self.parent.system.additional_data["combination_rule"] = combination_rule
        self.parent.system.additional_data["gen_pairs"] = gen_pairs == "yes"
        self.parent.system.additional_data["fudgeLJ"] = fudgeLJ
        self.parent.system.additional_data["fudgeQQ"] = fudgeQQ
        tokens.assert_no_more_than(5)

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
