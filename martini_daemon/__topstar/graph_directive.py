from typing import Any

import difflib
from ..__parser import (
    register_directive,
    GromacsTopFile,
    Directive,
    TokenList,
    TokenParseException,
)
from ..__core import System
from .graph import Graph, GraphAtomType

keyword_to_GraphAtomType = {
    "atom": GraphAtomType.NORMAL,
    "atom?": GraphAtomType.OPT,
    "atom!": GraphAtomType.NOT,
}


@register_directive
class GraphDirective(Directive):
    def __init__(self, parent: GromacsTopFile, path, line_num) -> None:
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
                self.graph.atoms.append((part_id, name_pat, type_pat, type_))
            case "equivalent":
                parts = set(tokens.unwrap(i, "word") for i in range(1, len(tokens)))
                self.graph.equivalents.append(parts)
            case _:
                if keyword in self.filters:
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

    def finish(self):
        self.graph.finish_init()
        self.graphs[self.graph.name] = self.graph

    @classmethod
    def is_mandatory(cls):
        return False

    @classmethod
    def is_unique(cls):
        return False

    @classmethod
    def is_valid_parent(cls, parent: Any) -> bool:
        return isinstance(parent, GromacsTopFile)

    @classmethod
    def get_name(cls) -> str:
        return "graph"

    aliases = {"frag", "fragment"}
