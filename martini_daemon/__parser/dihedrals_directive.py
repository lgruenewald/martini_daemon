from collections.abc import Callable

from ..__core import Force
from .gromacs_top_file import register_directive
from .interaction_directive import InteractionDirective


@register_directive
class DihedralsDirective(InteractionDirective):
    @classmethod
    def get_number_members(cls) -> int:
        return 4

    __type_data: dict[int, tuple[str, list[str]]] = {}

    @classmethod
    def register_type(
        cls, class_: type[Force], type_: int, name: str, args: list[str]
    ) -> None:
        assert type_ not in cls.__type_data
        cls.__type_data[type_] = (name, args)

    @classmethod
    def get_type(cls, type_int: int) -> str | None:
        got = cls.__type_data.get(type_int)
        return got[0] if got is not None else None

    @classmethod
    def get_type_args(cls, type_int: int) -> list[str]:
        got = cls.__type_data.get(type_int)
        assert got is not None
        return got[1]

    @classmethod
    def get_name(cls) -> str:
        return "dihedrals"

    @classmethod
    def get_number_params(cls, type_: int) -> tuple[int, int]:
        _, args = cls.__type_data[type_]
        return len(args), len(args)


def register_dihedral_type(
    type_: int, args: list[str]
) -> Callable[[type[Force]], type[Force]]:
    def inner(class_: type[Force]) -> type[Force]:
        name = class_.get_name()
        DihedralsDirective.register_type(class_, type_, name, args)
        return class_

    return inner
