import math
import re

from .token import Token


class TokenParseException(Exception):
    """Exception when .unwrap() fails"""

    def __init__(self, token: Token, message: str):
        self.token = token
        self.message = message


class TokenList:
    """A list of tokens.

    Parser calls the line method of :doc:`Directive</autoapi/martini_daemon/Directive>` with this type as the argument.
    Directive implementations should call the unwrap method with the
    """

    __int_pat = re.compile("^[-+]?[0-9]+$")
    __float_pat = re.compile("^[-+]?[0-9]+(\\.[0-9]*)?([eE][-+]?[0-9]+)?$")
    __word_pat = re.compile("^[a-zA-Z0-9_.]+$")
    __pattern_pat = re.compile("^[a-zA-Z0-9_?!*{}]+$")
    __pair_pat = re.compile(r"^[0-9]+:[a-zA-Z0-9_.]+$")
    __string_pat = re.compile(r'^"[^"]*"|<[^>]*>$')

    __DEFAULT = object()

    def __init__(self, line: str, tokens: list[Token], defines: dict[str, str]):
        """Note that tokenization has to performed first before constructing TokenList. This is done by the Parser.
        This is because the Parser needs to process line continuations first.

        :param line: Line that was tokenized.
        :param tokens: List of tokens obtained.
        """
        self.__line = line
        self.__tokens = tokens
        self.__defines = defines

    def get_line(self) -> str:
        return self.__line

    def __len__(self) -> int:
        return len(self.__tokens)

    def __getitem__(self, index: int) -> Token:
        return self.__tokens[index]

    def __setitem__(self, index: int, token: Token) -> None:
        self.__tokens[index] = token

    def assert_no_more_than(self, count: int) -> None:
        if len(self.__tokens) > count:
            raise TokenParseException(
                self.__tokens[count],
                f"Received unexpected additional tokens. Maximum number of tokens on this line is {count}.",
            )

    def assert_at_least(self, count: int) -> None:
        if len(self.__tokens) < count:
            raise TokenParseException(
                self.__tokens[count],
                f"Minimum number of tokens on this line is {count}.",
            )

    def unwrap(
        self,
        index: int,
        type_filter: str,
        default=__DEFAULT,
        error_msg: str | None = None,
    ):
        """Given a TokenList try to index it and convert to a usable value
        based on type_filter. If the index is out of range, a default
        value can be specified in place. Will also perform preprocessor #define
        replacements.

        :param index: Index of token in line.
        :param type_filter: What type to attempt to extract. See below.
        :param default: Default value to return if index is out of range.
        :param error_msg: The error message to print if the type filter is wrong.

        Possible filters:

        * int - returns int type, no processing.
        * float - returns float type, no processing.
        * positive - float, but raises an exception if 0 or smaller.
        * index - int, but subtracts 1, since GROMACS uses 1 based indexing, but OpenMM 0 based, exception if 0 or smaller.
        * degree - degree to radian and makes sure it's in the range -pi to pi, by adding or subtracting 2 pi.
        * word - a string with only alphanumerics and no whitespace in it.
        * pattern - fnmatch pattern, converts {} to [], as [] is used in fnmatch, but would conflict with directives.
        * pair - a tuple of an index (which reactant) and a word (which graph atom), separated by a colon.
        * string - "" or <> enclosed string. Strips the enclosing double quotes or angle parentheses.
        * raw - return the string content of the token, no validation. Useful to perform preprocessor #defines.
        """
        if len(self.__tokens) <= index:
            if default == self.__DEFAULT:
                # hack so "None" can also be used as a default value.
                # if it's still __DEFAULT it means no value was specified by the user.
                raise TokenParseException(
                    self.__tokens[-1],
                    f"Not enough tokens, expected token at index {index}.",
                )
            return default

        tok = self.__tokens[index]
        content = tok.line[tok.start : tok.end]

        while (got := self.__defines.get(content)) is not None:
            content = got

        match type_filter:
            case "int":
                if self.__int_pat.match(content):
                    return int(content)
            case "float":
                if self.__float_pat.match(content):
                    return float(content)
            case "positive":
                if self.__float_pat.match(content):
                    if float(content) <= 0.0:
                        raise TokenParseException(
                            tok,
                            error_msg or "Expected a positive non-zero real number.",
                        )
                    return float(content)
            case "index":
                if self.__int_pat.match(content):
                    if int(content) <= 0:
                        raise TokenParseException(
                            tok,
                            error_msg
                            or "Expected index, got an integer 0 or smaller."
                            "Note: indexing in .itp/.top files is usually 1 based.",
                        )
                    return int(content) - 1
            case "degree":
                if self.__float_pat.match(content):
                    angle = float(content) * math.pi / 180.0
                    while angle < math.pi:
                        angle += math.tau
                    while angle > math.pi:
                        angle -= math.tau
                    return angle
            case "word":
                if self.__word_pat.match(content):
                    return content
            case "pattern":
                if self.__pattern_pat.match(content):
                    return content.replace("{", "[").replace("}", "]")
            case "pair":
                # specialized index:word construct for [reaction] stuff
                if self.__pair_pat.match(content):
                    items = content.split(":")
                    return int(items[0]) - 1, items[1]
            case "string":
                if self.__string_pat.match(content):
                    return content[1:-1]
            case "raw":
                return content
            case _:
                raise TokenParseException(
                    tok, f"Filter {type_filter} couldn't be understood."
                )

        raise TokenParseException(
            tok, error_msg or f"Expected token of type {type_filter}."
        )
