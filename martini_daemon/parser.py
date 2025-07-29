# gromacs-style .ini format parser (used for .top/.itp files)

import re
from dataclasses import dataclass
from enum import Enum
import os
import math
from typing import Callable
import traceback


@dataclass
class Token:
    content: str
    path: str
    line: int
    start: int
    end: int


int_pat = re.compile("[-+]?[0-9]+")
float_pat = re.compile("[-+]?[0-9]+(\\.[0-9]*)?([eE][-+]?[0-9]+)?")
word_pat = re.compile("[a-zA-Z0-9_.]+")
pattern_pat = re.compile("[a-zA-Z0-9_?!*{}]+")
pair_pat = re.compile(r"[0-9]+:[a-zA-Z0-9_.]+")


class TokenParseError(Exception):
    def __init__(self, token: Token, message: str):
        self.token = token
        self.message = message


class ParseError(Exception):
    def __init__(self, message: str):
        self.message = message


def unwrap(tokens, index, type_filter, default="default placeholder"):
    """
    Given a list of tokens, try to index it and convert to a useable value
    based on type_filter. If the index would be out of range, a default
    value can be specified in place.

    Possible filters:
    int - returns int type, no processing
    float - returns float type, no processing
    positive - float, but raises an exception if 0 or smaller
    index - int, but subtracts 1, exception if 0 or smaller
    degree - degree to radian and makes sure it's in the range 0 to 2pi
    word - a string with only alphanumerics and no whitespace in it
    pattern - fnmatch pattern, converts {} to [] so {} can be used for sets
    pair - a tuple of an index (which reactant) and a word (which graph atom)
    """
    if len(tokens) <= index:
        if default == "default placeholder":
            # hack so "None" can also be used as a default value
            raise TokenParseError(
                tokens[-1],
                f"Not enough tokens, expected token at index {index}."
            )
        else:
            return default

    tok = tokens[index].content
    match type_filter:
        case "int":
            if int_pat.match(tok):
                return int(tok)
        case "float":
            if float_pat.match(tok):
                return float(tok)
        case "positive":
            if float_pat.match(tok):
                if float(tok) <= 0.:
                    raise TokenParseError(
                        tokens[index],
                        "Expected a positive non-zero real number."
                    )
                return float(tok)
        case "index":
            if int_pat.match(tok):
                if int(tok) <= 0:
                    raise TokenParseError(
                        tokens[index],
                        "Expected index, got an integer 0 or smaller."
                        "Note: indexing in .itp/.top files is usually 1 based."
                    )
                return int(tok) - 1
        case "degree":
            if float_pat.match(tok):
                angle = float(tok) * math.pi / 180.0
                while angle < 0.:
                    angle += math.tau
                while angle > math.tau:
                    angle -= math.tau
                return angle
        case "word":
            if word_pat.match(tok):
                return tok
        case "pattern":
            if pattern_pat.match(tok):
                return tok.replace("{", "[").replace("}", "]")
        case "pair":
            # specialized index:word construct for [reaction] stuff
            if pair_pat.match(tok):
                items = tok.split(":")
                return (int(items[0])-1, items[1])

    raise TokenParseError(
        tokens[index],
        f"Expected token of type {type_filter}."
     )


class Parser():
    """
    A parser for gromacs topology file-like config files with limited C
    style preprocessing, that should be sufficient for most input files.
    You can add your own levels (directives) as callbacks and this class will
    give you processed tokens for every line.

    You can parse the tokens using unwrap() from this module.

    All exceptions inside callbacks and otherwise will be caught and re-raised,
    so that the specific line and file they occured on can be printed to
    stdout, and the traceback surpressed to reduce clutter.

    Usage:
    p = Parser()
    p.add_level(name, handler)
    ok = p.parse(path_to_top, defines={})
    if not ok:
        # there was an error
    """

    def __init__(self):
        self._levels: dict[str, Callable] = {}
        self._start: dict[str, Callable] = {}
        self._end: dict[str, Callable] = {}
        self._mandatory: set[str] = set()
        self._unique: set[str] = set()

    def parse(self, path, include_dir, defines={}):
        """The main interface for using a TopParser class
        """
        self._linenum: int = 0
        self._path: str = path
        self._defines: dict[str, str] = defines
        self._include_dir: str | None = include_dir
        self._included: set[str] = set()
        self._past_directives: set[str] = set()
        self._current_level = None
        try:
            self._parse(path)
            for dir in self._mandatory:
                if dir not in self._past_directives:
                    raise ParseError(
                        f"Mandatory directive {dir} not found."
                    )
            return True
        except TokenParseError as pe:
            print("Parsing error:")
            print(pe.message)
            print()
            self.error_message_location(
                pe.token.path, pe.token.line,
                pe.token.start, pe.token.end
            )
            print()
        except ParseError as pe:
            print("Parsing error:")
            print(pe.message)
            print()
            self.error_message_location(self._path, self._linenum)
            print()
        except Exception:
            traceback.print_exc()
            print()
            print("Exception occured while parsing at:")
            self.error_message_location(self._path, self._linenum)
        return False

    def error_message_location(self, path, linenum, start=None, end=None):
        print(f"In file {path} at line {linenum+1}.")
        with open(path, "r") as f:
            line = f.read().split("\n")[linenum]
        if len(line) > 0:
            if (start is None or end is None):
                print(line)
            else:
                print(
                    f"{line[0:start]}"
                    f"\033[1;33m{line[start:end]}\033[0m"
                    f"{line[end:]}"
                )

    def _tokenize(self, line):
        """Splits a line up into a list of tokens. Similar to separating by
        whitespace, but more intelligent, and aware of things like ; comments.
        Also returns the type of tokens, and auto converts to int/float.
        """

        line = re.sub(r";.*", "", line).strip()
        toks = [
            Token(
                match.group(), self._path, self._linenum,
                match.start(), match.end()
            )
            for match in re.finditer(r"\S+", line)
        ]
        for i in range(len(toks)):
            while toks[i].content in self._defines:
                toks[i].content = self._defines[toks[i].content]
        return toks, line

    def _parse(self, path):
        """
        Parses path, adding new data linespython add lines to stack trace
        or directives to the accumulated list of directives so far.
        """

        self._included.add(path)
        oldpath = self._path
        self._path = path
        IfstackElem = Enum("IfstackElem", ["DoBranch", "SkipBranch",
                                           "SkippedIf", "Root"])
        ifstack = [IfstackElem.Root]  # ifstack is per file
        with open(path, "r") as fhandle:
            cumulative = ""
            for i, line in enumerate(fhandle):
                self._linenum = i
                # handle ignoring line endings
                if len(line) > 0 and line[-1] == "\\":
                    cumulative += line[:-1]
                    continue
                elif len(cumulative) > 0:
                    line = cumulative + line
                    cumulative = ""
                tokens, line = self._tokenize(line)

                # empty lines ignored
                # ignored only after checking for \ -- means that empty lines
                # also need \ to keep continuing one long line
                if len(tokens) == 0:
                    continue

                ifstack_top = ifstack[-1]
                # handle if/else logic before other things
                if tokens[0].content[0] == "#":
                    match tokens[0].content:
                        case "#ifdef":
                            if len(tokens) != 2:
                                raise TokenParseError(
                                    tokens[0], "#ifdef takes one argument."
                                )
                            if ifstack_top in {IfstackElem.DoBranch,
                                               IfstackElem.Root}:
                                if (
                                    float_pat.match(tokens[1].content)
                                    and float(tokens[1].content) != 0.0
                                ):
                                    ifstack.append(IfstackElem.DoBranch)
                                else:
                                    ifstack.append(IfstackElem.SkipBranch)
                            else:
                                ifstack.append(IfstackElem.SkippedIf)
                            continue

                        case "#ifndef":
                            if len(tokens) != 2:
                                raise TokenParseError(
                                    tokens[0], "#ifndef takes one argument."
                                )
                            if ifstack_top in {IfstackElem.DoBranch,
                                               IfstackElem.Root}:
                                if (
                                    float_pat.match(tokens[1].content)
                                    and float(tokens[1].content) != 0.0
                                ):
                                    ifstack.append(IfstackElem.SkipBranch)
                                else:
                                    ifstack.append(IfstackElem.DoBranch)
                            else:
                                ifstack.append(IfstackElem.SkippedIf)
                            continue
                        case "#else":
                            if len(tokens) != 1:
                                raise TokenParseError(
                                    tokens[0], "#else takes no argument."
                                )
                            if ifstack_top == IfstackElem.DoBranch:
                                ifstack[-1] = IfstackElem.SkipBranch
                            elif ifstack_top == IfstackElem.SkipBranch:
                                ifstack[-1] = IfstackElem.DoBranch
                            elif ifstack_top == IfstackElem.Root:
                                raise TokenParseError(
                                    tokens[0], "#else unmatched."
                                )
                            continue
                        case "#endif":
                            if len(tokens) != 1:
                                raise TokenParseError(
                                    tokens[0], "#endif takes no argument."
                                )
                            if ifstack_top == IfstackElem.Root:
                                raise TokenParseError(
                                    tokens[0], "#endif unmatched."
                                )
                            ifstack.pop()
                            continue
                        case "#end":
                            raise TokenParseError(
                                tokens[0], "Please use #endif."
                            )
                        # must list all other macro words here
                        # so that it doesn't error

                        # pass these here because only evaluate the if/else
                        # things here...
                        case "#include":
                            pass
                        case "#define":
                            pass
                        case "#undef":
                            pass
                        case _:
                            # error at unknown macros
                            raise TokenParseError(
                                tokens[0], "Unknown preprocessor directive."
                            )

                # everything below this only happens if the #ifdef/#else
                # says it should happen
                if ifstack_top in {IfstackElem.SkipBranch,
                                   IfstackElem.SkippedIf}:
                    continue

                if tokens[0].content[0] == "[":
                    # directives
                    if line[-1] != "]":
                        raise TokenParseError(
                            tokens[-1], "Invalid directive, not closed by ']'."
                        )
                    if tokens[0].content == "[":
                        tok_index = 1
                    else:
                        tok_index = 0
                    end_hook = self._end.get(self._current_level)
                    end_hook and end_hook()
                    self._current_level = line.strip("[] \t")
                    if (
                        self._current_level in self._unique
                        and self._current_level in self._past_directives
                    ):
                        raise ParseError(
                            tokens[tok_index],
                            "Unique directive present more than once."
                        )
                    self._past_directives.add(self._current_level)
                    start_hook = self._start.get(self._current_level)
                    start_hook and start_hook()
                    if self._levels.get(self._current_level) is None:
                        raise TokenParseError(
                            tokens[tok_index],
                            f"Unknown directive {self._current_level}."
                        )
                elif tokens[0].content[0] == "#":
                    match tokens[0].content:
                        case "#include":
                            name = line.replace("#include", "").strip()
                            if name[0] == '"' and name[-1] == '"':
                                name = name.strip('"')
                            elif name[0] == "<" and name[-1] == ">":
                                name = name.strip("<>")
                            else:
                                raise TokenParseError(
                                    Token(
                                        "", self._path, self._linenum,
                                        len("#include "), len(line)
                                    ),
                                    "#include argument should be inside"
                                    " quotation marks or <>."
                                )
                            search_dirs = [os.path.dirname(path),
                                           self._include_dir]
                            found = False
                            for dir in search_dirs:
                                newpath = os.path.join(dir, name)
                                if os.path.isfile(newpath):
                                    if newpath in self._included:
                                        raise TokenParseError(
                                            Token(
                                                "", self._path, self._linenum,
                                                len("#include "), len(line)
                                            ),
                                            f"Double inclusion of {newpath}."
                                        )
                                    self._parse(newpath)
                                    found = True
                                    break
                            if not found:
                                raise TokenParseError(
                                    Token(
                                        "", self._path, self._linenum,
                                        len("#include "), len(line)
                                    ),
                                    f"File not found: {name}."
                                )
                        case "#define":
                            if len(tokens) not in {2, 3}:
                                raise TokenParseError(
                                    tokens[0],
                                    "#define takes one or two arguments."
                                    " Preprocessor macros are not supported."
                                )
                            key = tokens[1].content
                            val = "1"
                            if len(tokens) == 3:
                                if tokens[2].content[0] in '"<':
                                    raise TokenParseError(
                                        tokens[2],
                                        "#define with strings/paths as tokens"
                                        " is not supported."
                                    )
                                val = tokens[2].content
                            self._defines[key] = val
                        case "#undef":
                            if len(tokens) != 2:
                                raise TokenParseError(
                                    tokens[0],
                                    "#undef takes one argument."
                                )
                            key = tokens[1].content
                            if self._defines.get(key):
                                self._defines.pop(key)
                        case _:
                            raise TokenParseError(
                                tokens[0],
                                f"Unknown preprocessor directive {tokens[0].content}."
                            )
                else:
                    # data lines
                    if self._current_level == "":
                        raise TokenParseError(
                            tokens[0],
                            "Data line encountered outside of any directive."
                        )
                    handler = self._levels.get(self._current_level)
                    if handler is None:
                        raise TokenParseError(
                            tokens[0],
                            "Data line encountered in unknown directive."
                        )
                    handler(tokens)

        if len(ifstack) > 1:
            raise ParseError(
                "Unmatched #ifdef or #ifndef."
                " #ifdef/#ifndef crossing file boundaries are not supported."
            )
        self._path = oldpath

    def add_level(
            self, name, handler, start=None, end=None,
            mandatory=False, unique=False
        ):
        """
        Add a new level to this Parser.

        name is the string contained within brackets / directive name.

        handler will be called on every line (with tokens, a dataclass
        containing the string content, file path and line number).
        In principle, the tokens are whitespace separated strings of a single
        line, but this handles backslash escaping of newlines and simple
        #defines (no preprocessor macros, only "variables", which can be
        recursive).

        start is a callback called every time at the start of a directive,
        end is a callback called every time at the end of a directive.
        """
        self._levels[name] = handler
        if start is not None:
            self._start[name] = start
        if end is not None:
            self._end[name] = end
        if mandatory:
            self._mandatory.add(name)
        if unique:
            self._unique.add(name)
