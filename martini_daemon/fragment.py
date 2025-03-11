from typing import Optional
from collections import OrderedDict


class Fragment():
    name: str
    particles: OrderedDict[str, int]
    opt: OrderedDict[str, Optional[int]]
    frag_id: int

    def __init__(self, name, id):
        self.name = name
        self.particles = OrderedDict()
        self.opt = OrderedDict()
        self.frag_id = id
        self.graph = None

    def get_atom(self, id: int | str) -> tuple[bool, Optional[int]]:
        if type(id) is int:
            id = f"{id+1}"
        if id in self.particles:
            return (True, self.particles.get(id))
        elif id in self.opt:
            # None is a valid value
            return (True, self.opt.get(id))
        else:
            return (False, None)

    # len and index_atom together are used in the iterator Fragments
    # this is meant to represent the numbered version of the fragment,
    # so it does not include the optional atoms for constant length
    def len(self) -> int:
        return len(self.particles)

    def index_atom(self, i: int) -> int:
        # TODO this could be slow, memoize if slow
        return list(self.particles.values())[i]


class Fragments():
    # a view over a fragment, particles or multiple fragments that
    # can be indexed
    # TODO: use this for all types of list[Fragment]
    # TODO: remove int based indexing for consistency
    # TODO try to remove parts for simplicity/consistency
    frags: list[Fragment]
    parts: list[int]

    def __init__(self, frags: list[Fragment], parts: list[int]):
        self.frags = frags
        self.parts = parts

    def from_parts(parts: list[int]):
        return Fragments([], parts)

    def from_frag(frag: Fragment):
        return Fragments([frag], [])

    def from_frags(frags: list[Fragment]):
        return Fragments(frags, [])

    def get(self, id: int | str | tuple[int, str | int]
            ) -> tuple[bool, Optional[int]]:
        if type(id) in {str, int}:
            if len(self.frags) == 0 and len(self.parts) > 0:
                if type(id) is str:
                    raise NotImplementedError
                return True, self.parts[id]
            elif len(self.parts) == 0 and len(self.frags) == 1:
                return self.frags[0].get_atom(id)
            else:
                raise ValueError("int/str index on a Fragments with wrong len")
        assert len(self.parts) == 0
        index, part_name = id
        return self.frags[index].get_atom(part_name)

    def particles(self):
        """Returns an iterator over all non-optional particles in child fragments."""
        class FragmentsIterator():
            def __init__(self, fragments: Fragments):
                self.frags = fragments

            def __iter__(self):
                self.cfrag = 0
                self.cindex = -1
                return self

            def __next__(self) -> tuple[int, str]:
                if len(self.frags.parts) > 0 and len(self.frags.frags) == 0:
                    self.cindex += 1
                    if self.cindex >= len(self.frags.parts):
                        raise StopIteration
                    return self.frags.parts[self.cindex]
                elif len(self.frags.frags) > 0 and len(self.frags.parts) == 0:
                    self.cindex += 1
                    if self.cindex >= self.frags.frags[self.cfrag].len():
                        self.cfrag += 1
                        self.cindex = 0
                    if self.cfrag >= len(self.frags.frags):
                        raise StopIteration
                    return \
                        self.frags.frags[self.cfrag].index_atom(self.cindex)

        return FragmentsIterator(self)
