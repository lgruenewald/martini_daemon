#!/usr/bin/env python3
# python class that groups together all the template information
# TODO update this

from parser import TopParser
from dataclasses import dataclass
import sys


@dataclass
class FuncAtom:
    rid: int
    role: str


@dataclass
class FuncMapAtom:
    rid: int
    id: int
    type: str
    atomname: str


class FuncMap:
    moltype: str
    atoms: list[FuncMapAtom]

    def __init__(self):
        pass

    def __str__(self):
        return "FuncMap for " + self.moltype + "{\n" + \
            "\n".join([str(i) for i in self.atoms]) + "\n}"

    def parse(self, dirs):
        dir = dirs.prevdir()
        if len(dir.data) != 1 or len(dir.data[0]) != 1 or \
                dir.data[0][0].type != "word":
            dirs.error("1 line 1 word token expected in [fmap] in [fgroup]")
            return False
        self.moltype = dir.data[0][0].content
        if not dirs.match("atoms"):
            dirs.error("[atoms] expected inside [fmap]")
            return False
        self.atoms = []
        atoms = dirs.prevdir()
        for line in atoms.data:
            if len(line) != 4:
                dirs.error("4 tokens expected per line in [atoms] in [fmap]")
                return False
            rid = line[0]
            id = line[1]
            type = line[2]
            name = line[3]
            match rid.type, id.type, type.type, name.type:
                case "int", "int", "word", "word":
                    self.atoms.append(FuncMapAtom(rid.content,
                                                  id.content,
                                                  type.content,
                                                  name.content))
                case _:
                    dirs.error("int, int, word, word expected in [atoms]")
                    return False
        return True


class FuncGroup:
    name: str
    atoms: list[FuncAtom]
    fmaps: list[FuncMap]

    def __init__(self):
        pass

    def __str__(self):
        return "FuncGroup " + self.name + "\natoms:{\n" + \
            "\n".join([str(i) for i in self.atoms]) + "\n}\n" +\
            "\n".join([str(i) for i in self.fmaps])

    def parse(self, dirs):
        dir = dirs.prevdir()
        if len(dir.data) != 1 or len(dir.data[0]) != 1 or \
                dir.data[0][0].type != "word":
            dirs.error("1 line of 1 word token expected in [fgroup] directive")
            return False
        self.name = dir.data[0][0].content
        if not dirs.match("atoms"):
            dirs.error("[atoms] expected inside [fgroup]")
            return False
        atoms = dirs.prevdir()
        self.atoms = []
        for line in atoms.data:
            if len(line) != 2:
                dirs.error("2 tokens expected in [atoms] in [fgroup]")
                return False
            rid = line[0]
            role = line[1]
            match rid.type, role.type:
                case "int", "word":
                    self.atoms.append(FuncAtom(rid.content, role.content))
                case _:
                    dirs.error("Int, word expected")
                    return False
        self.fmaps = []
        while dirs.match("fmap"):
            fm = FuncMap()
            if not fm.parse(dirs):
                return False
            self.fmaps.append(fm)
        return True


@dataclass
class ReactionMapAtom:
    reactant_id: str
    reactant_atom: int
    product_id: str
    product_atom: int


class Reaction:
    name: str
    r1: str
    r2: str
    p1: str
    distance_max: float
#    constraints TODO
    mapping: list[ReactionMapAtom]

    def __init__(self):
        pass

    def __str__(self):
        return "Reaction " + self.name + \
               "\nr1: " + self.r1 + " r2: " + self.r2 + " p1: " + self.p1 + \
               " distance_max: " + str(self.distance_max) + "\n" + \
               "\n".join([str(i) for i in self.mapping])

    def parse(self, dirs):
        dir = dirs.prevdir()
        # key val pairs in [reaction]
        for line in dir.data:
            if len(line) != 3:
                dirs.error("3 token long lines expected", line[0])
                return False
            key = line[0]
            eq = line[1]
            value = line[2]
            if eq.content != "=":
                dirs.error("second token must be '='", line[1])
                return False
            match value.type, key.content:
                case "word", "name":
                    self.name = value.content
                case "word", "r1":
                    self.r1 = value.content
                case "word", "r2":
                    self.r2 = value.content
                case "word", "p1":
                    self.p1 = value.content
                case "float" | "int", "distance_max":
                    self.distance_max = value.content
                case _:
                    dirs.error("Invalid combo of key name and val type", key)
                    return False
        # TODO constraints
        if not dirs.match("constraints"):
            dirs.error("[constraints] expected in [reaction]")
            return False
        # mapping
        if not dirs.match("mapping"):
            dirs.error("[mapping] expected in [reaction] after constraints")
            return False
        mapdir = dirs.prevdir()
        self.mapping = []
        for line in mapdir.data:
            if len(line) != 6:
                dirs.error("6 token long lines expected", line[0])
                return False
            rname = line[0]
            colon1 = line[1]
            rid = line[2]
            pname = line[3]
            colon2 = line[4]
            pid = line[5]
            if colon1.content != ":" or colon2.content != ":":
                dirs.error("Tokens 2 and 5 need to be ':'", colon1)
                return False
            match rname.type, rid.type, pname.type, pid.type:
                case "word", "int", "word", "int":
                    self.mapping.append(ReactionMapAtom(rname.content,
                                                        rid.content,
                                                        pname.content,
                                                        pid.content))
                case _:
                    dirs.error("word, int, word, int expected", rname)
                    return False
        return True


class Template:
    funcs: list[FuncGroup]
    reactions: list[Reaction]

    def __init__(self):
        self.funcs = []
        self.reactions = []

    def __str__(self):
        return "\n".join([str(i) for i in self.funcs] +
                         [str(i) for i in self.reactions])

    def parse(self, path, defines):
        p = TopParser()
        ok, dirs = p.parse(path, defines)
        if not ok:
            return False
        while dirs.peek() != "":
            if dirs.match("fgroup"):
                fg = FuncGroup()
                if not fg.parse(dirs):
                    return False
                self.funcs.append(fg)
            elif dirs.match("reaction"):
                rx = Reaction()
                if not rx.parse(dirs):
                    return False
                self.reactions.append(rx)
            # TODO reaction product
            else:
                tok = dirs.peek_token()
                dirs.error("Unknown directive", tok)
                return False
        return True


if __name__ == "__main__":
    argv = sys.argv
    argc = len(argv)
    if argc != 2:
        print("Usage: ./template.py <file>")
        quit(1)
    path = argv[1]

    t = Template()
    ok = t.parse(path, {})
    print("OK?", ok)
    if ok:
        print(t)
