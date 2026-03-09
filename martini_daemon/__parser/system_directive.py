from typing import Any

from .directive import Directive
from .gromacs_top_file import GromacsTopFile, register_directive
from .token_list import TokenList

@register_directive
class SystemDirective(Directive):
    def __init__(self, parent: GromacsTopFile, path: str, line_num: int) -> None:
        self.parent = parent
        self.path = path
        self.line_num = line_num

    def where(self) -> tuple[str, int]:
        return self.path, self.line_num

    def line(self, tokens: TokenList) -> None:
        self.parent.set("title", tokens.get_line())

    def finish(self):
        pass

    @classmethod
    def is_mandatory(cls):
        return True

    @classmethod
    def is_unique(cls):
        return True

    @classmethod
    def is_valid_parent(cls, parent: Any) -> bool:
        return type(parent) == GromacsTopFile

    @classmethod
    def get_name(cls) -> str:
        return "system"