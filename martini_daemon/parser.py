#!/usr/bin/env python3
# Martini Topology + fragments and reaction templates parser

import re
import sys
from dataclasses import dataclass
from enum import Enum
from collections import OrderedDict
import os
import math


@dataclass
class Token:
    content: 'typing.Any'
    type: str
    path: str
    line: int


def unwrap(tokens, index, type_filter, default="default placeholder"):
    if len(tokens) <= index:
        if default == "default placeholder":
            # hack so "None" can also be used as a default value
            raise ValueError(f"Not enough tokens, expect token index {index}")
        else:
            return default
    set_filter_pass = type(type_filter) is set and tokens[index].type in type_filter
    if tokens[index].type != type_filter and not set_filter_pass:
        if tokens[index].type == "int" and type_filter == "float":
            # the only implicit conversion we do is int -> float
            return float(tokens[index].content)
        elif tokens[index].type == "int" and type_filter == "index":
            # the other magic done here is subtracting 1 from indices
            return tokens[index].content - 1
        elif ((tokens[index].type == "int" or tokens[index].type == "float")
              and type_filter == "degree"):
            # magic 3: convert degrees to radians
            return tokens[index].content * math.pi / 180.0
        raise ValueError(f"Token {tokens[index].content}: "
                         f"expected {type_filter} "
                         f"but received {tokens[index].type} instead")
    else:
        return tokens[index].content


class TopParser:
    """A parser for gromacs topology file-like config files

    Usage:
    p = TopParser()
    p.add_level(name, handler)
    ok, dirs = p.parse(path_to_top, defines={})
    if ok:
        # do stuff with dirs
    """

    # line number for reporting errors
    _linenum: int
    # current file being parsed
    _path: str
    # was there an error at any point?
    _haderror: bool
    # which files were #included
    _included: dict[str, bool]
    # which DEFINES were defined
    _defines: dict
    # current level similar to how openmm parses top files
    _current_level = None
    # TODO potential improvement: validate the nesting of levels to see if the
    # structure makes sense. It wouldn't change behavior on correct top files
    # but would enhance error messages on incorrect top files
    _include_dir: str

    # list of functions to call with data lines in each level
    _levels: dict
    # list of directives we already complained about
    _complained_directives: dict

    def __init__(self):
        self._complained_directives = {}
        self._levels = {}

    def parse(self, path, include_dir, defines={}):
        """The main interface for using a TopParser class
        """
        self._linenum = 0
        self._path = path
        self._haderror = False
        self._defines = defines
        self._include_dir = include_dir
        self._complained_directives = {}
        self._included = {}
        self._parse(path)
        return not self._haderror

    def error(self, message):
        print("\033[1;91mParser error\033[0m")
        print(f" {message}")
        print(f" In file {self._path} at line {self._linenum + 1}")
        self._haderror = True

    def _tokenize(self, line):
        """Splits a line up into a list of tokens. Similar to separating by
        whitespace, but more intelligent, and aware of things like ; comments.
        Also returns the type of tokens, and auto converts to int/float.
        """

        # token patterns and their types
        patterns = OrderedDict()
        patterns["int"] = re.compile("[-+]?[0-9]+")
        patterns["float"] = re.compile("[-+]?[0-9]+(\\.[0-9]*)?"
                                       "([eE][-+]?[0-9]+)?")
        patterns["key"] = re.compile("[a-zA-Z0-9_]+:")
        patterns["word"] = re.compile("[a-zA-Z0-9_.]+")
        patterns["pattern"] = re.compile("[a-zA-Z0-9_?!*{}]+")
        patterns["macro"] = re.compile("#[a-zA-Z0-9_]+")
        patterns["string"] = re.compile('"[^"]*"')
        patterns["bracket_string"] = re.compile("<[^>]*>")
        patterns["key"] = re.compile(r"[a-zA-Z0-9_]+\s*:")
        patterns["symbol"] = re.compile(r"[\[\]:='.?!+-]")

        cur = 0  # current position in the line
        tokens = []  # list of tokens built up so far
        key = None  # TODO this is so horrible

        while cur < len(line):
            ch = line[cur]
            if ch in " \t\r\n":
                cur += 1  # ignore whitespace
            elif ch == ";":
                break  # ignore comments (until end of line)
            else:
                # find the right pattern
                # the one that matches the most characters
                # is considered the right pattern
                # if multiple match the equal length, the first to match
                # (according to the order in patterns) will be the one
                longest_type = None
                longest_span = 0
                for type, pat in patterns.items():
                    match = pat.match(line, cur)
                    if match is None:
                        continue
                    start, end = match.span()
                    # only longest matches than previously possible should
                    # overwrite matches that came earlier
                    if end - start > longest_span:
                        longest_span = end - start
                        longest_type = type
                # no matching token
                if longest_type is None:
                    if self._current_level != "system":
                        # allow anything in [system] since that's used
                        # for basically a comment
                        self.error("Unexpected character " + ch)
                    cur += 1
                    continue
                content = line[cur:cur + longest_span]
                # processing of content -> convert to numbers, remove quotes
                if longest_type == "float":
                    content = float(content)
                elif longest_type == "int":
                    content = int(content)
                elif longest_type == "string":
                    content = content.strip('"')
                elif longest_type == "bracket_string":
                    content = content.strip("<>")
                elif longest_type == "pattern":
                    content = content.replace("{", "[").replace("}", "]")
                elif longest_type == "key":
                    if key is not None:
                        raise ValueError("Double key")
                    key = content.strip(": \t")
                    cur += longest_span
                    continue
                # #define value replacements (only for word tokens)
                while longest_type == "word":
                    val = self._defines.get(content)
                    # todo more elegant preservation of types
                    if val is None:
                        break
                    match val:
                        case int():
                            longest_type = "int"
                        case float():
                            longest_type = "float"
                        case str():
                            longest_type = "word"
                        case _:
                            self.error("Unknown DEFINE type, can't replace")
                    content = val
                if key is not None:
                    content = (key, content)
                    longest_type = "pair"
                    key = None
                # make token
                tok = Token(content, longest_type,
                            self._path, self._linenum)
                tokens.append(tok)
                cur += longest_span

        return tokens

    def _parse(self, path):
        """Parses path, adding new data linespython add lines to stack trace or directives to the accumulated
        list of directives so far.
        """

        # stuff for #includes and the #ifdef stack
        if path in self._included:
            self.error("Double inclusion of " + path)
            return
        self._included[path] = True
        oldpath = self._path
        self._path = path
        IfstackElem = Enum("IfstackElem", ["DoBranch", "SkipBranch",
                                           "SkippedIf", "Root"])
        ifstack = [IfstackElem.Root]  # ifstack is per file
        with open(path, "r") as fhandle:
            cumulative = []
            for i, line in enumerate(fhandle):
                self._linenum = i
                tokens = self._tokenize(line)
                # handle ignoring line endings
                if len(tokens) > 0 and tokens[-1].content == "\\":
                    cumulative += tokens[:-1]
                    continue
                elif len(cumulative) > 0:
                    tokens = cumulative + tokens
                    cumulative = []

                # empty lines ignored
                # ignored only after checking for \ -- means that empty lines
                # also need \ to keep continuing one long line
                if len(tokens) == 0:
                    continue

                ifstack_top = ifstack[-1]
                # handle if/else logic before other things
                if tokens[0].type == "macro":
                    match tokens[0].content:
                        case "#ifdef":
                            if len(tokens) != 2:
                                self.error("#ifdef takes one argument")
                                continue
                            if ifstack_top in {IfstackElem.DoBranch,
                                               IfstackElem.Root}:
                                if tokens[1].type in {"int", "float"} and \
                                        tokens[1].content > 0:
                                    ifstack.append(IfstackElem.DoBranch)
                                else:
                                    ifstack.append(IfstackElem.SkipBranch)
                            else:
                                ifstack.append(IfstackElem.SkippedIf)
                            continue

                        case "#ifndef":
                            if len(tokens) != 2:
                                self.error("#ifndef takes one argument")
                                continue
                            if ifstack_top in {IfstackElem.DoBranch,
                                               IfstackElem.Root}:
                                if tokens[1].type in {"int", "float"} and \
                                        tokens[1].content > 0:
                                    ifstack.append(IfstackElem.SkipBranch)
                                else:
                                    ifstack.append(IfstackElem.DoBranch)
                            else:
                                ifstack.append(IfstackElem.SkippedIf)
                            continue
                        case "#else":
                            if len(tokens) != 1:
                                self.error("#else takes no argument")
                            if ifstack_top == IfstackElem.DoBranch:
                                ifstack[-1] = IfstackElem.SkipBranch
                            elif ifstack_top == IfstackElem.SkipBranch:
                                ifstack[-1] = IfstackElem.DoBranch
                            elif ifstack_top == IfstackElem.Root:
                                self.error("#else unmatched")
                            continue
                        case "#endif":
                            if len(tokens) != 1:
                                self.error("#endif takes no argument")
                            if ifstack_top == IfstackElem.Root:
                                self.error("#endif unmatched")
                                continue
                            ifstack.pop()
                            continue
                        case "#end":
                            self.error("Please use #endif")
                            continue
                        # must list all other macro words here
                        # so that it doesn't error
                        case "#include":
                            pass
                        case "#define":
                            pass
                        case "#undef":
                            pass
                        case _:
                            # error at unknown macros
                            self.error("Unknown " + tokens[0].content)
                if ifstack_top in {IfstackElem.SkipBranch,
                                   IfstackElem.SkippedIf}:
                    continue

                if tokens[0].content == "[":
                    # directives
                    if len(tokens) != 3:
                        self.error("Invalid directive: wrong len(tokens)")
                        continue
                    if tokens[2].content != "]":
                        self.error("Invalid directive: no ]")
                    self._current_level = tokens[1].content
                    if self._levels.get(self._current_level) is None:
                        if self._complained_directives.get(self._current_level):
                            self._haderror = True
                        else:
                            self.error(f"Unknown directive: {self._current_level}")
                            self._complained_directives[self._current_level] = True
                elif tokens[0].type == "macro":
                    match tokens[0].content:
                        case "#include":
                            if len(tokens) != 2:
                                self.error("#include should have 1 argument")
                            if tokens[1].type not in {"string",
                                                      "bracket_string"}:
                                self.error("#include argument should be inside"
                                           "quotations or <>")
                            name = tokens[1].content
                            search_dirs = [os.path.dirname(path),
                                           self._include_dir]
                            found = False
                            for dir in search_dirs:
                                newpath = os.path.join(dir, name)
                                if os.path.isfile(newpath):
                                    self._parse(newpath)
                                    found = True
                                    break
                            if not found:
                                print(tokens)
                                self.error(f"File not found: {name}")
                        case "#define":
                            if not (len(tokens) in {2, 3}) or \
                                    tokens[1].type != "word":
                                self.error("#define takes one word argument")
                                self.error("or 1 word arg and 1 any arg")
                                continue
                            key = tokens[1].content
                            val = 1
                            if len(tokens) == 3:
                                if tokens[2].type not in \
                                        {"word", "int", "float"}:
                                    self.error("Can only #define numbers "
                                               "and words")
                                val = tokens[2].content
                            self._defines[key] = val
                        case "#undef":
                            if len(tokens) != 2 or tokens[1].type != "word":
                                self.error("#undef takes one word argument")
                                continue
                            key = tokens[1].content
                            if self._defines.get(key):
                                self._defines.pop(key)
                        case _:
                            self.error("Unknown " + tokens[0].content)
                else:
                    # data lines
                    if self._current_level == "":
                        self.error("Data line outside of directives")
                        continue
                    handler = self._levels.get(self._current_level)
                    if handler is None:
                        if self._complained_directives.get(self._current_level):
                            self._haderror = True
                        else:
                            self.error(f"Data line in unknown directive "
                                       f"{self._current_level}")
                            self._complained_directives[self._current_level] = True
                        continue
                    try:
                        handler(tokens)
                    except Exception as e:
                        print("\033[1;91mCallback error\033[0m")
                        print(f"In file {path} at line {self._linenum}")
                        raise e
        if len(ifstack) > 1:
            self.error("Unmatched #ifdef or #ifndef")
        self._path = oldpath

    def add_level(self, name, handler):
        """Add a new level to this TopParser
        """
        self._levels[name] = handler


if __name__ == "__main__":
    argv = sys.argv
    argc = len(argv)
    if argc != 2:
        print("Usage: ./top_parser.py <file>")
        quit(1)
    path = argv[1]
    p = TopParser()
    ok, dir = p.parse(path)
    print("OK?", ok)
    if ok:
        print(dir)
