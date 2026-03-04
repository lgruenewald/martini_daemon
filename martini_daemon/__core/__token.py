from dataclasses import dataclass

@dataclass
class Token:
    content: str
    path: str
    line: int
    start: int
    end: int


