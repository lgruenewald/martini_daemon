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
    content: str
#    type: str
    path: str
    line: int


int_pat = re.compile("[-+]?[0-9]+")
float_pat = re.compile("[-+]?[0-9]+(\\.[0-9]*)?([eE][-+]?[0-9]+)?")
word_pat = re.compile("[a-zA-Z0-9_.]+")
pattern_pat = re.compile("[a-zA-Z0-9_?!*{}]+")


def unwrap(tokens, index, type_filter, default="default placeholder"):
    if len(tokens) <= index:
        if default == "default placeholder":
            # hack so "None" can also be used as a default value
            raise ValueError(f"Not enough tokens, expect token index {index}")
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
        case "index":
            if int_pat.match(tok):
                return int(tok) - 1
        case "degree":
            if float_pat.match(tok):
                return float(tok) * math.pi / 180.0
        case "word":
            if word_pat.match(tok):
                return tok
        case "pattern":
            if pattern_pat.match(tok):
                return tok.replace("{", "[").replace("}", "]")

    raise ValueError(f"Token {tok}: "
                     f"expected {type_filter}.")


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
    _start: dict
    _end: dict
    # list of directives we already complained about
    _complained_directives: dict

    def __init__(self):
        self._complained_directives = {}
        self._levels = {}
        self._start = {}
        self._end = {}

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

        line = re.sub(r";.*", "", line).strip()
        toks = [Token(tok, self._path, self._linenum) for tok in line.split()]
        for i in range(len(toks)):
            while toks[i].content in self._defines:
                toks[i].content = self._defines[toks[i].content]
        return toks, line

    def _parse(self, path):
        """Parses path, adding new data linespython add lines to stack trace
        or directives to the accumulated list of directives so far.
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
                                self.error("#ifdef takes one argument")
                                continue
                            if ifstack_top in {IfstackElem.DoBranch,
                                               IfstackElem.Root}:
                                if float_pat.match(tokens[1].content) and float(tokens[1].content) != 0.0:
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
                                if float_pat.match(tokens[1].content) and float(tokens[1].content) != 0.0:
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
                            self.error("Unknown preprocessor directive "
                                       f"{tokens[0].content}")

                # everything below this only happens if the #ifdef/#else
                # says it should happen
                if ifstack_top in {IfstackElem.SkipBranch,
                                   IfstackElem.SkippedIf}:
                    continue

                if tokens[0].content[0] == "[":
                    # directives
                    if line[-1] != "]":
                        self.error("Invalid directive: no ]")
                    end_hook = self._end.get(self._current_level)
                    end_hook and end_hook()
                    self._current_level = line.strip("[] \t")
                    start_hook = self._start.get(self._current_level)
                    start_hook and start_hook()
                    if self._levels.get(self._current_level) is None:
                        if self._complained_directives.get(self._current_level):
                            self._haderror = True
                        else:
                            self.error(f"Unknown directive: {self._current_level}")
                            self._complained_directives[self._current_level] = True
                elif tokens[0].content[0] == "#":
                    match tokens[0].content:
                        case "#include":
                            name = line.replace("#include", "").strip()
                            if name[0] == '"' and name[-1] == '"':
                                name = name.strip('"')
                            elif name[0] == "<" and name[-1] == ">":
                                name = name.strip("<>")
                            else:
                                self.error("#include argument should be inside"
                                           "quotations or <>")
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
                                self.error(f"File not found: '{name}'")
                        case "#define":
                            if len(tokens) not in {2, 3}:
                                self.error("#define takes one or two arguments")
                                continue
                            key = tokens[1].content
                            val = "1.0"
                            if len(tokens) == 3:
                                if tokens[2][0] in '"<':
                                    self.error("Cannot #define strings/paths")
                                val = tokens[2].content
                            self._defines[key] = val
                        case "#undef":
                            if len(tokens) != 2:
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

    def add_level(self, name, handler, start=None, end=None):
        """Add a new level to this TopParser
        """
        self._levels[name] = handler
        if start is not None:
            self._start[name] = start
        if end is not None:
            self._end[name] = end
