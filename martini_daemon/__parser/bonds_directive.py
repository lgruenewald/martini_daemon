from .interaction_directive import InteractionDirective
from .gromacs_top_file import register_directive

@register_directive
class BondsDirective(InteractionDirective):
    @classmethod
    def get_number_members(cls) -> int:
        return 2

    __type_data: dict[int, tuple[str, list[str]]] = {}
    __is_exclusion: dict[str, bool] = {}

    @classmethod
    def register_type(cls, type_: int, name: str, args: list[str], is_excl: bool) -> None:
        cls.__type_data[type_] = (name, args)
        cls.__is_exclusion[name] = is_excl

    @classmethod
    def get_type(cls, type_int: int) -> str | None:
        got = cls.__type_data.get(type_int)
        return got or got[0]

    @classmethod
    def get_type_args(cls, type_int: int) -> list[str]:
        return cls.__type_data.get(type_int)[1]

    @classmethod
    def get_name(cls) -> str:
        return "bonds"

    @classmethod
    def is_exclusion(cls, type_: str) -> bool:
        return cls.__is_exclusion[type_]

# TODO import and add types to class_
def register_bond_type(class_, type_: int, args: list[str], is_excl: bool) -> None:
    name = class_.get_name()
    BondsDirective.register_type(type_, name, args, is_excl)
    return class_