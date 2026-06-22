from dataclasses import dataclass


@dataclass(frozen=True)
class Token:
    """
    A single whitespace separated token from parsing.

    :param content: The processed content of the token (before preprocessing, TokenList.unwrap() runs macro substitutions).
    :param line_num: The line number in the source file.
    :param path: The path to the source file.
    :param start: The start position in the line for this Token (before preprocessing).
    :param end: The end position in the line for this Token (before preprocessing).
    """

    content: str
    line_num: int
    path: str
    start: int
    end: int
