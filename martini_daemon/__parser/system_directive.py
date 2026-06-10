from .directive import Directive
from .gromacs_top_file import GromacsTopFile, register_directive
from .token_list import TokenList
from .parser import ParseException


@register_directive
class SystemDirective(Directive):
    def line(self, tokens: TokenList) -> None:
        assert isinstance(self.parent, GromacsTopFile)
        if title := self.parent.system.additional_data.get("title") is not None:
            raise ParseException(f"Double title definition, title was already defined to be {title}")
        self.parent.system.additional_data["title"] = tokens.get_line()

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
        return "system"
