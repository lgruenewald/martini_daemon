from typing import Any
from abc import ABCMeta, abstractmethod

from .directive import Directive
from .token_list import TokenList, TokenParseException
from .molecule_type_directive import MoleculeTypeDirective


class InteractionDirective(Directive, metaclass=ABCMeta):
    # default Directive boilerplate
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
        return isinstance(parent, MoleculeTypeDirective)

    # get_name should be implemented by child classes

    # customizable line parsing
    def line(self, tokens: TokenList) -> None:
        type_num, type_name = self.read_type(tokens)

        self.parent.interactions.append((
            type_name,
            self.read_members(tokens),
            self.read_params(tokens, type_num),
            self.is_exclusion(type_name),
        ))

    # reasonable defaults that still can be overridden for e.g. virtual_sitesn
    def read_members(self, tokens: TokenList) -> list[int]:
        return [
            self.parent.parse_index(tokens, i)
            # 0 to n
            for i in range(self.get_number_members())
        ]

    def read_params(self, tokens: TokenList, type_: int) -> list[float]:
        min_params, max_params = self.get_number_params(type_)
        tokens.assert_no_more_than(
            self.get_number_members() + 1 + max_params
        )
        tokens.assert_at_least(
            self.get_number_members() + 1 + min_params
        )
        params_start = self.get_number_members() + 1
        params = []
        for j, filter_ in enumerate(self.get_type_args(type_)):
            # ignore optional params
            if params_start + j >= len(tokens):
                break
            params.append(
                tokens.unwrap(
                    params_start + j,
                    filter_,
                )
            )

        return params

    def read_type(self, tokens: TokenList) -> tuple[int, str]:
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
            return type_num, type_


    # for most child classes, number of members leads to a good default impl of read_type, read_members, read_params
    @classmethod
    @abstractmethod
    def get_number_members(cls) -> int:
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def get_number_params(cls, type_: int) -> tuple[int, int]:
        """
        Return the minimum and maximum number of params.
        """
        raise NotImplementedError

    # metadata for molecule type
    @classmethod
    def is_exclusion(cls, type_: str) -> bool:
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



