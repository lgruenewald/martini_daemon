from dataclasses import dataclass


@dataclass(frozen=True)
class Token:
    """A single whitespace separated token from parsing.

    :param content: The processed content of the token.
    :param line: The whole line that contains the token. Before preprocessing, after removing comments and making lines whole across backslashes.
    :param line_num: The line number in the source file.
    :param path: The path to the source file.
    :param start: The start position in the line for this Token.
    :param end: The end position in the line for this Token.
    """

    content: str
    line: str
    line_num: int
    path: str
    start: int
    end: int
