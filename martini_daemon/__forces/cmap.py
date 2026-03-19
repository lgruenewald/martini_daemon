from typing import Any
import openmm as mm

from ..__core import BondedForce, register_available_force
from ..__parser import InteractionDirective, register_directive, Directive, GromacsTopFile, TokenList, TokenParseException

@register_directive
class CMAPTypeDirective(Directive):
    def line(self, tokens: TokenList) -> None:
        parts = tuple(
            self.parent.unwrap_atom_type(tokens, i)
            for i in range(5)
        )
        type_ = tokens.unwrap(5, "int")
        if type_ != 1:
            raise TokenParseException(
                tokens[5],
                f"Unsupported CMAP type {type_}"
            )
        size = tokens.unwrap(6, "int")
        if size != tokens.unwrap(7, "int"):
            raise TokenParseException(
                tokens[7],
                "Non-square CMAPs are not supported."
            )
        if size < 8:
            raise TokenParseException(
                tokens[6],
                "CMAPs of size below 8 are too small. This might result in large errors in interpolation (compared to GROMACS)."
            )

        params = [
            tokens.unwrap(8+i, "float") for i in range(size*size)
        ]

        # rearrangement as in
        # https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/gromacstopfile.py
        # lines 959-984
        # based on testing:
        # gromacs CMAPs start at -180,-180
        # rows are the second dihedral, going from -180 to +180
        # columns are the first dihedral, going from -180 to +180

        cmap = []
        midpoint = size // 2
        for n in range(size):
            column = (n + midpoint) % size
            for o in range(size):
                row = (o + midpoint) % size
                cmap.append(params[size * row + column])


        cmaps = self.parent.system.additional_data.get("cmap_maps")
        types = self.parent.system.additional_data.get("cmap_types")
        if types is None:
            types = {}
            self.parent.system.additional_data["cmap_types"] = types
        if cmaps is None:
            cmaps = []
            self.parent.system.additional_data["cmap_maps"] = cmaps

        types[parts] = len(cmaps)
        cmaps.append((size, cmap))

    def finish(self):
        pass

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
        return "cmaptypes"


@register_directive
class CMAPDirective(InteractionDirective):
    @classmethod
    def get_number_params(cls, type_: int) -> tuple[int, int]:
        return 0, 0

    @classmethod
    def get_number_members(cls) -> int:
        return 5

    @classmethod
    def get_type(cls, type_int: int) -> str | None:
        if type_int != 1:
            return None
        return "cmap"

    @classmethod
    def get_type_args(cls, type_int: int) -> list[str]:
        return []

    @classmethod
    def get_name(cls) -> str:
        return "cmap"

    @classmethod
    def is_exclusion(cls, type_: str) -> bool:
        return False

@register_available_force
class Cmap(BondedForce):

    def _add_to_force(self, members: list[int], params: list[float]) -> None:
        assert len(params) == 0
        types = tuple(
            self.system.get_type(member) for member in members
        )
        map = self.system.additional_data["cmap_types"].get(types)
        if map is None:
            raise ValueError(f"Unknown CMAP type for {types}.")

        i, j, k, l, m = members
        self.force.addTorsion(
            map,
            i, j, k, l,
            j, k, l, m
        )


    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        return params

    @staticmethod
    def uses_pbc() -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "cmap"

    def _set_force_obj(self):
        self.force = mm.CMAPTorsionForce()
        for index, (size, cmap) in enumerate(self.system.additional_data["cmap_maps"]):
            i = self.force.addMap(size, cmap)
            assert i == index