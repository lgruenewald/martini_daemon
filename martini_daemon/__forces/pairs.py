import openmm as mm

from ..__parser import Directive, register_directive, GromacsTopFile, TokenList, InteractionDirective, TokenParseException
from ..__core import BondedForce, register_available_force

@register_directive
class PairTypes(Directive):
    def line(self, tokens: TokenList) -> None:
        pair_types = self.parent.system.additional_data.get("pairtypes")
        if pair_types is None:
            pair_types = {}
            self.parent.system.additional_data["pairtypes"] = pair_types
        t1 = self.parent.unwrap_atom_type(tokens, 0)
        t2 = self.parent.unwrap_atom_type(tokens, 1)
        type_ = tokens.unwrap(2, "int")
        if type_ != 1:
            raise TokenParseException(
                tokens[2],
                f"Unsupported pair type {type_}."
            )
        sigma = tokens.unwrap(3, "float")
        epsilon = tokens.unwrap(4, "float")
        pair_types[(t1, t2)] = (sigma, epsilon)

    def finish(self):
        pass

    @classmethod
    def is_mandatory(cls):
        return False

    @classmethod
    def is_unique(cls):
        return False

    @classmethod
    def is_valid_parent(cls, parent: GromacsTopFile) -> bool:
        return isinstance(parent, GromacsTopFile)

    @classmethod
    def get_name(cls) -> str:
        return "pairtypes"

@register_directive
class PairsDirective(InteractionDirective):
    @classmethod
    def get_number_members(cls) -> int:
        return 2

    @classmethod
    def get_number_params(cls, type_: int) -> tuple[int, int]:
        return 0, 2

    @classmethod
    def get_type(cls, type_int: int) -> str | None:
        if type_int != 1:
            return None
        return "pair"

    @classmethod
    def get_type_args(cls, type_int: int) -> list[str]:
        return ["float", "float"]

    @classmethod
    def get_name(cls) -> str:
        return "pairs"

@register_available_force
class Pairs(BondedForce):
    def __init__(self, system):
        super().__init__(system)
        self.epsilon_r = system.additional_data["epsilon_r"]

    def _add_to_force(self, members: list[int], params: list[float]) -> None:
        self.force.addBond(*members, params)

    def _parse(self, members: list[int], params: list[float]) -> list[float]:
        # TODO what if charge/type changes during sim
        q1 = self.system.get_charge(members[0])
        q2 = self.system.get_charge(members[1])
        q_prod = q1 * q2

        if len(params) == 0:
            pair_types = self.system.additional_data.get("pairtypes")
            if pair_types is None:
                raise ValueError("No pair types were defined.")
            t1 = self.system.get_type(members[0])
            t2 = self.system.get_type(members[1])
            params = pair_types.get((t1, t2)) or pair_types.get((t2, t1))
            if params is None:
                raise ValueError(
                    f"Unknown pair type ({t1}, {t2})."
                    f" Valid types are: {pair_types.keys()}."
                )

        sigma, epsilon = params

        c6 = 4 * epsilon * (sigma ** 6)
        c12 = 4 * epsilon * (sigma ** 12)

        return [q_prod, c6, c12]


    @staticmethod
    def uses_pbc() -> bool:
        return True

    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def get_name(cls) -> str:
        return "pair"

    def _set_force_obj(self):
        self.force = mm.CustomBondForce(
            "LJ + ES;"
            "LJ = (C12 / r^12 - C6 / r^6);"
            "ES = f*qprod/epsilon_r/r;"
            f"epsilon_r = {self.epsilon_r};"
            "f = 138.935458;"
        )
        self.force.addPerBondParameter("qprod")
        self.force.addPerBondParameter("C6")
        self.force.addPerBondParameter("C12")

    def flag_atom_change(self, atom_id, change_charge) -> None:
        # TODO
        raise NotImplementedError
