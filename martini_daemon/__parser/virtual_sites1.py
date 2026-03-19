from .interaction_directive import InteractionDirective
from .gromacs_top_file import register_directive
from ..__parser import TokenList

@register_directive
class VirtualSites1(InteractionDirective):
    @classmethod
    def get_number_members(cls) -> int:
        return 2

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
        return "virtual_sites1"

    @classmethod
    def get_number_params(cls, type_: int) -> tuple[int, int]:
        _, args = cls.__type_data[type_]
        return len(args), len(args)

    def line(self, tokens: TokenList) -> None:
        # we don't want to trigger nrexcl processing with these exclusions
        super().line(tokens)
        self.parent.exclusions.add((
            self.parent.parse_index(tokens, 0),
            self.parent.parse_index(tokens, 1)
        ))

def register_vsite1_type(type_: int, args: list[str]):
    def inner(class_):
        name = class_.get_name()
        VirtualSites1.register_type(type_, name, args)
        return class_
    return inner
