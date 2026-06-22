import difflib

from ..__core import System
from ..__parser import (
    Directive,
    GromacsTopFile,
    TokenList,
    TokenParseException,
    register_directive,
)
from .graph import Graph, GraphAtomType

keyword_to_GraphAtomType = {
    "atom": GraphAtomType.NORMAL,
    "atom?": GraphAtomType.OPT,
    "atom!": GraphAtomType.NOT,
}


@register_directive
class GraphDirective(Directive):
    def __init__(self, parent: GromacsTopFile, path: str, line_num: int) -> None:
        super().__init__(parent, path, line_num)
        self.graph = Graph(None)
        self.system: System = parent.system

        if self.system.additional_data.get("graphs") is None:
            self.system.additional_data["graphs"] = {}
        self.graphs = self.system.additional_data["graphs"]
        self.filters = self.system.get_filters()

    def line(self, tokens: TokenList) -> None:
        keyword = tokens.unwrap(0, "raw")
        match keyword:
            case "name":
                if self.graph.name is not None:
                    raise TokenParseException(
                        tokens[0], f"Graph already has a name: {self.graph.name}."
                    )
                name = tokens.unwrap(1, "word")
                if self.graphs.get(name) is not None:
                    raise TokenParseException(
                        tokens[1], f"Redefinition of graph with name {name}."
                    )
                self.graph.name = name
            case "atom" | "atom?" | "atom!":
                part_id = tokens.unwrap(1, "word")
                name_pat = tokens.unwrap(2, "pattern")
                type_pat = tokens.unwrap(3, "pattern")
                type_ = keyword_to_GraphAtomType[keyword]
                if self.graph.atom_name_to_index.get(part_id) is not None:
                    raise TokenParseException(
                        tokens[1], f"Same graph has already a node called {part_id}."
                    )
                self.graph.add_atom(part_id, name_pat, type_pat, type_)
            case "equivalent":
                parts = set()
                for i in range(1, len(tokens)):
                    tok = tokens.unwrap(i, "word")
                    if self.graph.atom_name_to_index.get(tok) is None:
                        raise TokenParseException(
                            tokens[i],
                            f"Node {tok} was not yet defined. Define it before `equivalent`.",
                        )
                    if tok in parts:
                        raise TokenParseException(
                            tokens[i],
                            f"Node {tok} is specified twice on the same `equivalent` line.",
                        )
                    parts.add(tok)
                self.graph.equivalents.append(parts)
            case _:
                if keyword in self.filters:
                    parts = []
                    last_special = None
                    for i in range(1, len(tokens)):
                        tok = tokens.unwrap(i, "word")
                        idx = self.graph.atom_name_to_index.get(tok)
                        if idx is None:
                            raise TokenParseException(
                                tokens[i],
                                f"Interaction {keyword} references undefined node {tok}.",
                            )
                        if self.graph.atoms[idx][3] != GraphAtomType.NORMAL:
                            if last_special is not None:
                                raise TokenParseException(
                                    tokens[i],
                                    f"Interaction {keyword} references more than one "
                                    + f" optional or forbidden atoms ({last_special}, {tok}). "
                                    + "Forbidden/optional atoms connected to eachother are not allowed.",
                                )
                            last_special = tok
                        parts.append(tok)
                    parts = [tokens.unwrap(i, "word") for i in range(1, len(tokens))]
                    self.graph.interactions.append((keyword, parts))
                else:
                    possibilities = self.filters | {
                        "atom",
                        "atom?",
                        "atom!",
                        "name",
                        "equivalent",
                    }

                    close_matches = difflib.get_close_matches(keyword, possibilities, 3)
                    raise TokenParseException(
                        tokens[0],
                        f"Keyword {keyword} not recognized."
                        + (
                            f" Perhaps you meant one of: {', '.join(close_matches)}"
                            if len(close_matches) > 0
                            else ""
                        )
                        + f" Valid keywords are: {', '.join(possibilities)}.",
                    )

    def finish(self) -> None:
        self.graph.finish_init()
        self.graphs[self.graph.name] = self.graph

    @classmethod
    def is_mandatory(cls) -> bool:
        return False

    @classmethod
    def is_unique(cls) -> bool:
        return False

    @classmethod
    def is_valid_parent(cls, parent: Directive) -> bool:
        return isinstance(parent, GromacsTopFile)

    @classmethod
    def get_name(cls) -> str:
        return "graph"

    aliases = {"frag", "fragment"}
