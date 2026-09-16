from ..__core import System
from .directive import Directive
from .parser import Parser
from .token_list import TokenList, TokenParseException


class InvalidTopologyError(Exception):
    pass


class GromacsTopFile(Directive):
    # static fields
    __top_directives: list[type[Directive]] = []

    @classmethod
    def add_top_directive(cls, directive_: type[Directive]) -> None:
        """Add a directive to the global .top format parser."""
        cls.__top_directives.append(directive_)

    def __init__(
        self,
        system: System,
        path: str,
        include_dirs: list[str] | None = None,
        defines: dict[str, str] | None = None,
    ) -> None:
        """
        Parse a Gromacs Top File.

        This class is:
        * Metadata about the .top format, allowing for the construction of parsers and generators of said format.
        * A class you can instantiate given a .top file, which will perform the parsing into itself that you can then read out.

        :param path: Path to the .top file.
        :param system: __core.System.
        :param include_dirs: List of directories to be searched if #include fails to find a file in the current dir.
        :param defines: dict[str, str] of keys and values for token replacements by the limited C preprocessor impl.
        """
        if defines is None:
            defines = {}
        defines["DAEMON"] = ""
        # Parser's init should be called
        self.__parser = Parser(
            root=self, path=path, include_dirs=include_dirs, defines=defines
        )
        # Directive's __init__ is raise NotImplementedError
        for directive_ in self.__top_directives:
            self.__parser.add_directive(directive_)

        # results of parsing
        self.system = system

        ok = self.__parser.parse()
        if not ok:
            # error was already printed, but we re-raise here
            raise InvalidTopologyError

    def unwrap_atom_type(self, tokens: TokenList, index: int) -> str:
        type_ = tokens.unwrap(index, "word")
        if self.system.get_atom_type(type_) is None:
            raise TokenParseException(tokens[index], f"Unknown atom type {type_}.")
        return type_

    # === implementing Directive ===
    def line(self, tokens: TokenList) -> None:
        raise TokenParseException(
            tokens[0], "Data line encountered outside of any directive."
        )

    def finish(self) -> None:
        pass

    @classmethod
    def get_name(cls) -> str:
        return "<root GromacsTopFile>"

    # the other methods of Directive should not get called on the root, so it's fine
    @classmethod
    def is_mandatory(cls) -> bool:
        raise NotImplementedError()

    @classmethod
    def is_unique(cls) -> bool:
        raise NotImplementedError()

    @classmethod
    def is_valid_parent(cls, parent: Directive) -> bool:
        raise NotImplementedError()

    def where(self) -> tuple[str, int]:
        raise NotImplementedError()


def register_directive[T: type[Directive]](class_: T) -> T:
    GromacsTopFile.add_top_directive(class_)
    return class_
