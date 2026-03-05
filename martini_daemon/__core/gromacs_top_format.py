from typing import Type, Any
from .directive import Directive
from .parser import Parser
from .force import Force
from .vsite import VirtualSite
from .token_list import TokenList, TokenParseError
from .. import TokenList

class GromacsTopFormat(Directive):
    """
    This class is:
    * Metadata about the .top format, allowing for the construction of parsers and generators of said format.
    * A class you can instantiate given a .top file, which will perform the parsing into itself that you can then read out.
    """
    # static fields
    __top_directives: list[Type[Directive]] = []

    def __init__(
        self, path: str, include_dirs: list[str] | None = None,
            defines: dict[str, str] | None = None
    ):
        """
        :param path: Path to the .top file.
        :param include_dirs: List of directories to be searched if #include fails to find a file in the current dir.
        :param defines: dict[str, str] of keys and values for token replacements by the limited C preprocessor impl.
        """
        # Parser's init should be called
        self.__parser = Parser(
            root=self,
            path=path,
            include_dirs=include_dirs,
            defines=defines
        )
        # Directive's __init__ is raise NotImplementedError
        for directive_ in GromacsTopFormat.__top_directives:
            self.__parser.add_directive(directive_)

    # === implementing Directive ===
    def line(self, tokens: TokenList) -> None:
        raise TokenParseError(
            tokens[0],
            "Data line encountered outside of any directive."
        )

    # the other methods of Directive should not get called on the root, so it's fine
    def finish(self):
        raise NotImplementedError()

    @staticmethod
    def is_mandatory():
        raise NotImplementedError()

    @staticmethod
    def is_unique():
        raise NotImplementedError()

    @staticmethod
    def is_valid_parent(parent: Any) -> bool:
        raise NotImplementedError()

    @staticmethod
    def get_name() -> str:
        raise NotImplementedError()

    def where(self) -> tuple[str, int]:
        raise NotImplementedError()

    @staticmethod
    def add_top_directive(directive_: Type[Directive]) -> None:
        """
        Add a directive to the global .top format parser.
        """
        GromacsTopFormat.__top_directives.append(directive_)

def directive(class_: Type[Directive]) -> Type[Directive]:
    GromacsTopFormat.add_top_directive(class_)
    return class_

# these classes have decorators compatible with .force.Force
@directive
class __Exclusion(Directive):
    pass

def exclusion(class_: Type[Force]) -> Type[Force]:
    pass

@directive
class __Bond(Directive):
    pass

def bond(class_: Type[Force], type_: int) -> Type[Force]:
    pass

@directive
class __Angle(Directive):
    pass

def angle(class_: Type[Force], type_: int) -> Type[Force]:
    pass

@directive
class __Dihedral(Directive):
    pass

def dihedral(class_: Type[Force], type_: int) -> Type[Force]:
    pass

# these classes have decorators compatible with .vsite.VirtualSite
@directive
class __VirtualSite1(Directive):
    pass

def virtualsite1(class_: Type[VirtualSite], type_: int) -> Type[VirtualSite]:
    pass

@directive
class __VirtualSite2(Directive):
    pass

def virtualsite2(class_: Type[VirtualSite], type_: int) -> Type[VirtualSite]:
    pass

@directive
class __VirtualSite3(Directive):
    pass

def virtualsite3(class_: Type[VirtualSite], type_: int) -> Type[VirtualSite]:
    pass

@directive
class __VirtualSite4(Directive):
    pass

def virtualsite4(class_: Type[VirtualSite], type_: int) -> Type[VirtualSite]:
    pass

@directive
class __VirtualSiteN(Directive):
    pass

def virtualsiteN(class_: Type[VirtualSite], type_: int) -> Type[VirtualSite]:
    pass