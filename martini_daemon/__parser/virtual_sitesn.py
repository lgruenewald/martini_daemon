from .token_list import TokenList, TokenParseException
from .interaction_directive import InteractionDirective
from .gromacs_top_file import register_directive


@register_directive
class VirtualSitesN(InteractionDirective):
    @classmethod
    def get_number_members(cls) -> int:
        raise NotImplementedError

    @classmethod
    def get_number_params(cls, type_: int) -> tuple[int, int]:
        return len(cls.__type_data.get(type_)[1]), len(cls.__type_data.get(type_)[1])

    __type_data: dict[int, tuple[str, list[str]]] = {}

    @classmethod
    def register_type(cls, type_: int, name: str, args: list[str]) -> None:
        assert type_ not in cls.__type_data.keys()
        cls.__type_data[type_] = (name, args)

    @classmethod
    def get_type(cls, type_int: int) -> str | None:
        got = cls.__type_data.get(type_int)
        return got[0] if got is not None else None

    @classmethod
    def get_type_args(cls, type_int: int) -> list[str]:
        return cls.__type_data.get(type_int)[1]

    @classmethod
    def get_name(cls) -> str:
        return "virtual_sitesn"

    @classmethod
    def is_exclusion(cls, type_: str) -> bool:
        return False

    def read_members(self, tokens: TokenList) -> list[int]:
        n_params = self.get_number_params(tokens.unwrap(1, "int"))[0]
        start = 2
        end = len(tokens) - n_params
        step = n_params + 1
        return [self.parent.parse_index(tokens, 0)] + [
            self.parent.parse_index(tokens, i) for i in range(start, end, step)
        ]

    def read_params(self, tokens: TokenList, type_: int) -> list[float]:
        n_params = self.get_number_params(type_)[0]
        start = 2
        end = len(tokens) - n_params
        step = n_params + 1
        params = []
        for i in range(start, end, step):
            for j, filter_ in enumerate(self.get_type_args(type_)):
                params.append(
                    tokens.unwrap(
                        i + j + 1,
                        filter_,
                    )
                )

        return params

    def read_type(self, tokens: TokenList) -> tuple[int, str]:
        # vid type constructing...
        type_num = tokens.unwrap(1, "int")
        type_ = self.get_type(type_num)
        if type_ is None:
            raise TokenParseException(
                tokens[1], f"Invalid type {type_num} for directive {self.get_name()}"
            )
        else:
            return type_num, type_

    def line(self, tokens: TokenList) -> None:
        # we don't want to trigger nrexcl processing with these exclusions
        super().line(tokens)
        if len((members := self.read_members(tokens))) == 2:
            self.parent.molecule_type.exclusions.add(
                (
                    members[0],
                    members[1],
                )
            )


def register_vsiten_type(type_: int, args: list[str]):
    """
    args are per constructing atom
    """

    def inner(class_):
        name = class_.get_name()
        VirtualSitesN.register_type(type_, name, args)
        return class_

    return inner
