from __future__ import annotations

from abc import ABC, abstractmethod

from .token_list import TokenList


class Directive(ABC):
    """
    Base class for directives.

    You should not instantiate this directly.

    Directives should implement all of the methods.
    """

    def __init__(self, parent: Directive, path: str, line_num: int) -> None:
        """
        Start parsing a directive.

         Called each time a directive of this type occurs in the source file.

        :param parent: The instance of the parent directive or root object. Directives mutate this instance during parsing.
        :param path: Path to where it occurs.
        :param line_num: Line number where it occurs.
        """
        self.parent = parent
        self.path = path
        self.line_num = line_num

    def where(self) -> tuple[str, int]:
        """Return the path and line number passed in the constructor for error messages."""
        return self.path, self.line_num

    @abstractmethod
    def line(self, tokens: TokenList) -> None:
        """
        Parse a non-empty source line inside the directive.

        :param tokens: Preprocessed and tokenized source line.
        """
        raise NotImplementedError

    @abstractmethod
    def finish(self) -> None:
        """
        End directive parsing.

        This happens when a new directive is encountered, for which this directive is not a valid parent for.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def is_mandatory(cls) -> bool:
        """
        Return whether this directive is mandatory.

        Only valid for directives at root.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def is_unique(cls) -> bool:
        """
        Return whether this directive is unique.

        Only valid for directives at root.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def is_valid_parent(cls, parent: Directive) -> bool:
        """
        Return whether the provided object instance is of a valid type for parent of this directive.

        This is used to guarantee that certain directives are "inside" other directives,
        such as [bonds] being inside [moleculetype].

        Directives at root should "return True" for the root type.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def get_name(cls) -> str:
        """Return the name for the directive used in []."""
        raise NotImplementedError

    aliases = set()
