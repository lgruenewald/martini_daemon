#!/usr/bin/env python3
# generic topology-like file parser

import re
import sys
from dataclasses import dataclass


@dataclass
class Token:
    content: 'typing.Any'
    type: str
    path: str
    line: int


@dataclass
class Directive:
    name: str
    token: Token
    data: list[list[Token]]
    lines: list[str]


@dataclass
class DirectiveList:
    dirs: list[Directive]
    current: int = 0

    def prevdir(self):
        return self.dirs[self.current-1]

    def peek(self):
        if self.current < len(self.dirs):
            return self.dirs[self.current].name
        else:
            return ""

    def peek_token(self):
        if self.current < len(self.dirs):
            return self.dirs[self.current].token
        else:
            return Token("", "", "", 0)

    def match(self, name):
        if (self.peek() == name):
            self.advance()
            return True
        else:
            return False

    def advance(self):
        name = self.peek()
        if self.peek() != "":
            self.current += 1
        return name

    def error(self, message, token: Token):
        print(message)
        print("In file ", token.path, " at line ", token.line + 1)


class TopParser:

    # line number for reporting errors
    _linenum = 0
    # current file being parsed
    _path = ""
    # was there an error at any point?
    _haderror = False
    # which files were #included
    _included = {}
    # which DEFINES were defined
    _defines = {}

    _directives = []

    def __init__(self):
        pass

    def _directive(self, directive, token):
        self._directives.append(Directive(directive, token, [], []))

    def _data(self, tokens, line):
        # tokens - list of tokens, already tokenized and converted to number
        # line - raw string of the whole line, as in the file
        # cindex - how many data lines were there before in this directive
        if len(self._directives) == 0:
            self.error("Data line outside of directives")
            return
        self._directives[-1].data.append(tokens)
        self._directives[-1].lines.append(line.strip())

    def error(self, message):
        print("Parse error: ", message)
        print("In file ", self._path, "at line ", self._linenum + 1)
        self._haderror = True

    def _tokenize(self, line):
        start = 0
        cur = 0
        tokens = []

        def advance():
            nonlocal cur
            cur += 1
            if cur <= len(line):
                return line[cur-1]
            else:
                return "\0"

        def peek():
            if cur < len(line):
                return line[cur]
            else:
                return "\0"

        def make_token(type):
            tok = Token(line[start:cur], type, self._path, self._linenum)
            tokens.append(tok)

        def make_token_re(regex, type):
            nonlocal start, cur
            match = regex.search(line, start)
            if match is None:
                self.error("Internal error: bad regex")
                return
            start, cur = match.span()
            make_token(type)

        number = re.compile("[-+]?[0-9]+(\\.[0-9]*)?([eE][-+]?[0-9]+)?")
        word = re.compile("[#a-zA-Z_][a-zA-Z0-9_]*")

        while cur < len(line):
            ch = advance()
            start = cur-1
            if ch in [" ", "\t", "\r", "\n"]:
                pass
            elif ch == "\"":
                while cur < len(line) and (ch := advance()) != "\"":
                    pass
                make_token("string")
                tokens[-1].content = tokens[-1].content.strip("\"")
            elif ch.isdigit() or \
                    (ch in "+-") and peek().isdigit():
                make_token_re(number, "int")
                t = tokens[-1].content
                if "." in t or "e" in t or "E" in t:
                    tokens[-1].content = float(t)
                    tokens[-1].type = "float"
                else:
                    tokens[-1].content = int(t)
            elif ch in ["_", "#"] or ch.isalpha():
                make_token_re(word, "macro" if ch == "#" else "word")
            elif ch == ";":
                break
            else:
                make_token("symbol")

        return tokens

    def _parse(self, path):
        if path in self._included:
            self.error("Double inclusion of " + path)
            return
        self._included[path] = True
        oldpath = self._path
        self._path = path
        with open(path, "r") as fhandle:
            cumulative = []
            for i, line in enumerate(fhandle):
                self._linenum = i
                tokens = self._tokenize(line)
                if len(tokens) == 0:
                    continue

                # handle ignoring line endings
                if tokens[-1].content == "\\":
                    cumulative += tokens[:-1]
                    continue
                elif len(cumulative) > 0:
                    tokens = cumulative + tokens
                    cumulative = []

                # directives
                if tokens[0].content == "[":
                    if len(tokens) != 3:
                        self.error("Invalid directive: wrong len(tokens)")
                    if tokens[2].content != "]":
                        self.error("Invalid directive: no ]")
                    self._directive(tokens[1].content, tokens[1])
                elif tokens[0].type == "macro":
                    match tokens[0].content:
                        case "#include":
                            if len(tokens) != 2:
                                self.error("#include should have 1 argument")
                            path = tokens[1].content
                            self._parse(path)
                        case _:
                            self.error("Unknown " + tokens[0].content)
                # data
                else:
                    self._data(tokens, line)
        self._path = oldpath

    def parse(self, path, defines={}):
        self._linenum = 0
        self._path = path
        self._defines = defines
        self._parse(path)
        return not self._haderror, DirectiveList(self._directives)


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
