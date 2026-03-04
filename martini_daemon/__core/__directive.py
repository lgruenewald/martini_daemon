from __future__ import annotations
from .__simulation import Simulation

class Directive:
    def __init__(self, sim: Simulation, parent: None | Directive) -> None:
        raise NotImplementedError

    def line(self) -> None:
        raise NotImplementedError

    def __del__(self):
        raise NotImplementedError

    @staticmethod
    def is_mandatory():
        raise NotImplementedError

    @staticmethod
    def is_unique():
        raise NotImplementedError