from collections.abc import Callable

from ..__core import Force
from ..__rust import Token
from .gromacs_top_file import register_directive
from .interaction_directive import InteractionDirective
from .molecule_type_directive import MoleculeTypeDirective


@register_directive
class PairsDirective(InteractionDirective):
    @classmethod
    def get_number_members(cls) -> int:
        return 2

    __type_data: dict[int, tuple[str, list[str]]] = {}

    @classmethod
    def register_type(cls, type_: int, name: str, args: list[str]) -> None:
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
        return "pairs"

    def get_types(self) -> dict[tuple, list[list[Token]]]:
        assert isinstance(self.parent, MoleculeTypeDirective)
        return self.parent.system.pair_types

    @classmethod
    def get_number_params(cls, type_: int) -> tuple[int, int]:
        _, args = cls.__type_data[type_]
        return len(args), len(args)

    def fallback_params(self, type_: int, types: list[str]) -> list[float] | None:
        if type_ != 1:
            return None
        assert isinstance(self.parent, MoleculeTypeDirective)
        if not self.parent.system.additional_data.get("gen_pairs"):
            return None
        fudgeLJ = self.parent.system.additional_data.get("fudgeLJ") or 1.0
        assert len(types) == 2
        # nonbond params, or applies combination rule if needed
        sigma, epsilon = self.parent.system.get_sigma_epsilon(*types)
        return [sigma, epsilon * fudgeLJ]


def register_pair_type(
    type_: int, args: list[str]
) -> Callable[[type[Force]], type[Force]]:
    def inner(class_: type[Force]) -> type[Force]:
        name = class_.get_name()
        PairsDirective.register_type(type_, name, args)
        return class_

    return inner
