from typing import Any
from abc import ABCMeta, abstractmethod

from .directive import Directive
from .token_list import TokenList, TokenParseException
from .molecule_type import MoleculeType


class InteractionDirective(Directive, metaclass=ABCMeta):
    # default Directive boilerplate
    def __init__(self, parent: MoleculeType, path: str, line_num: int) -> None:
        self.parent = parent
        self.path = path
        self.line_num = line_num

    def where(self) -> tuple[str, int]:
        return self.path, self.line_num

    def finish(self):
        pass

    @classmethod
    def is_mandatory(cls):
        return False

    @classmethod
    def is_unique(cls):
        return False

    @classmethod
    def is_valid_parent(cls, parent: Any) -> bool:
        return issubclass(type(parent), Directive)

    # get_name should be implemented by child classes

    # customizable line parsing
    def line(self, tokens: TokenList) -> None:
        type_ = self.read_type(tokens)

        self.parent.interactions.append((
            type_,
            self.read_members(tokens),
            self.read_params(tokens),
            self.is_exclusion(type_),
            self.is_virtual_site(type_),
            self.is_constraint(type_)
        ))

    # reasonable defaults that still can be overridden for e.g. virtual_sitesn
    def read_members(self, tokens: TokenList) -> list[int]:
        return [
            tokens.unwrap(
                i, "index"
            )
            # 0 to n
            for i in range(self.get_number_members())
        ]

    def read_params(self, tokens: TokenList) -> list[float]:
        return [
            tokens.unwrap(
                j, "float"
            )
            # skip 0 to n and n+1 (type)
            for j in range(self.get_number_members()+1, len(tokens))
        ]

    def read_type(self, tokens: TokenList) -> str:
        type_num = tokens.unwrap(
            self.get_number_members(),
            "int",
        )
        type_ = self.get_type(type_num)
        if type_ is None:
            raise TokenParseException(
                tokens[self.get_number_members()],
                f"Invalid type {type_num} for directive {self.get_name()}"
            )
        else:
            return type_


    # for most child classes, number of members leads to a good default impl of read_type, read_members, read_params
    @classmethod
    @abstractmethod
    def get_number_members(cls) -> int:
        raise NotImplementedError

    # metadata for molecule type
    @classmethod
    def is_exclusion(cls, type_: str) -> bool:
        return False

    @classmethod
    def is_virtual_site(cls, type_: str) -> bool:
        return False

    @classmethod
    def is_constraint(cls, type_: str) -> bool:
        return False


    # each class should have its own way of registering new types, if they'd like to provide it
    @classmethod
    @abstractmethod
    def get_type(cls, type_int: int) -> str | None:
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def get_type_args(cls, type_int: int) -> list[str]:
        raise NotImplementedError



