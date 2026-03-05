from typing import Any

from .gromacs_top_file import directive
from .directive import  Directive
from .token_list import TokenList


@directive
class Molecule(Directive):
    def __init__(self, parent: Any, path: str, line_num: int) -> None:
        pass

    def where(self) -> tuple[str, int]:
        pass

    def line(self, tokens: TokenList) -> None:
        pass

    def finish(self):
        pass

    @staticmethod
    def is_mandatory():
        pass

    @staticmethod
    def is_unique():
        pass

    @staticmethod
    def is_valid_parent(parent: Any) -> bool:
        pass

    @staticmethod
    def get_name() -> str:
        return "moleculetype"
# TODO add a "process_nrel_excl" here instead of adding exclusions together with bonds