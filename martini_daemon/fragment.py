from typing import Optional
from collections import OrderedDict
from enum import Enum
from .graph import GraphFragment


class GotAtom(Enum):
    NotFound = 1
    MissingOptional = 2
    Found = 3


class Fragment():
    name: str
    particles: OrderedDict[str, int]
    opt: OrderedDict[str, Optional[int]]
    frag_id: int
    graph: Optional[GraphFragment]

    def __init__(self, name, id):
        self.name = name
        self.particles = OrderedDict()
        self.opt = OrderedDict()
        self.frag_id = id
        self.graph = None

    def get_atom(self, id: str) -> tuple[GotAtom, Optional[int]]:
        if type(id) is int:
            id = f"{id+1}"
        if id in self.particles:
            return (GotAtom.Found, self.particles.get(id))
        elif id in self.opt:
            # None is a valid value, signaling a missing optional
            opt_got = self.opt.get(id)
            got_atom = opt_got is None and GotAtom.MissingOptional or GotAtom.Found
            return (got_atom, opt_got)
        else:
            return (GotAtom.NotFound, None)

    # len and index_atom together are used in the iterator Fragments
    # this is meant to represent the numbered version of the fragment,
    # so it does not include the optional atoms for constant length
    def len(self) -> int:
        return len(self.particles)

    def index_atom(self, i: int) -> int:
        # TODO this could be slow, memoize if slow
        return list(self.particles.values())[i]


def index_pair(
                frags: list[Fragment] | list[int], pair: int | tuple[int, str]
              ) -> tuple[GotAtom, Optional[int]]:
    if type(pair) is int and type(frags) is list:
        if pair >= 0 and pair < len(frags):
            return (GotAtom.Found, frags[pair])
        else:
            return (GotAtom.NotFound, None)
    elif type(pair) is tuple and type(frags) is list:
        return frags[pair[0]].get_atom(pair[1])
    else:
        raise ValueError(f"bad type for frags {type(frags)} pair {type(pair)}")
