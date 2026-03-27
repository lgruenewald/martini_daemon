from typing import Any
import difflib
from math import pi, cos

from ..__parser import Directive, register_directive, GromacsTopFile, TokenList, ParseException, TokenParseException, MoleculeTypeDirective
from ..__rust import DetectionTemplate
from .modification_template import ModificationTemplate
from .graph import GraphAtomType


@register_directive
class ReactionDirective(MoleculeTypeDirective):

    def __init__(self, parent, path, line_num):
        super().__init__(parent, path, line_num)
        self.d_template = DetectionTemplate()
        self.molecule_type: ModificationTemplate = ModificationTemplate()
        self.molecule_type.nrexcl = 1
        if parent.system.additional_data.get("detection_templates") is None:
            parent.system.additional_data["detection_templates"] = []
        parent.system.additional_data["detection_templates"].append(self.d_template)
        self.system = parent.system


    def line(self, tokens: TokenList) -> None:
        name = tokens.unwrap(0, "word")
        self.d_template.name = name
        self.molecule_type.name = name
        tokens.assert_no_more_than(1)
        if self.system.molecule_types.get(name) is not None:
            raise TokenParseException(
                tokens[0],
                f"Reaction name {name} conflicts with existing reaction or molecule name."
            )


    def finish(self):
        super().finish()
        try:
            # called during parsing
            self.d_template.complete()
        except Exception as e:
            raise ParseException(str(e))

    @classmethod
    def is_mandatory(cls):
        return False

    @classmethod
    def is_unique(cls):
        return False

    @classmethod
    def is_valid_parent(cls, parent: Any) -> bool:
        return isinstance(parent, GromacsTopFile)

    aliases = {"rx"}
    @classmethod
    def get_name(cls) -> str:
        return "reaction"

    def parse_pair(self, tokens: TokenList, index: int) -> tuple[int, int]:
        if len(self.d_template.reactants) == 0:
            raise ParseException(
                "[reactants] must come before being able to index atoms within a reaction."
            )
        r_index, a_name = tokens.unwrap(index, "pair")
        if r_index < 0 or r_index >= len(self.d_template.reactants):
            raise TokenParseException(
                tokens[index],
                f"Reactant index {r_index} is out of range (1 to {len(self.d_template.reactants)})"
            )
        graph = self.system.additional_data["graphs"][self.d_template.reactants[r_index]]
        a_index = graph.atom_name_to_index.get(a_name)
        if a_index is None:
            raise TokenParseException(
                tokens[index],
                f"Reaction atom {a_name} not found in graph {graph.name}."
            )
        if graph.atoms[a_index][3] == GraphAtomType.NOT:
            raise TokenParseException(
                tokens[index],
                f"Reaction atom {a_name} is a forbidden atom, so it is impossible to reference it."
            )
        return r_index, a_index


    # returns a flattened view, since MoleculeType does things in a flattened way
    def parse_index(self, tokens: TokenList, index: int) -> int:
        r_index, a_index = self.parse_pair(tokens, index)
        flattened_index = a_index
        for i in range(r_index):
            # all reactants that come before
            graph = self.system.additional_data["graphs"][self.d_template.reactants[i]]
            flattened_index += len(graph.atoms)
        return flattened_index

@register_directive
class ReactantsDirective(Directive):
    def line(self, tokens: TokenList) -> None:
        if len(self.parent.d_template.reactants) > 0:
            raise TokenParseException(
                tokens[0],
                "Only one line containing all reactants per reaction."
            )
        reactants = [
            tokens.unwrap(i, "word") for i in range(len(tokens))
        ]
        for i, r in enumerate(reactants):
            if (
                self.parent.system.additional_data.get("graphs") is None
                or self.parent.system.additional_data["graphs"].get(r) is None
            ):
                raise TokenParseException(
                    tokens[i],
                    f"Undefined reactant graph name {r}. Please put the required [graph] directive first."
                )
        self.parent.molecule_type.reactants = reactants
        self.parent.d_template.reactants = reactants

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
        return isinstance(parent, ReactionDirective)

    aliases = {"reactant"}
    @classmethod
    def get_name(cls) -> str:
        return "reactants"

@register_directive
class ConditionsDirective(Directive):

    @staticmethod
    def __parse_angle_conditions(tokens: TokenList, start: int, wrap: bool) -> list[tuple[int, int]]:
        """
        If wrap is false, the "angle space" is 0 to pi
        If wrap is true, the "angle space" is -pi to pi and is considered periodic
        """
        # 0 to 20 or 50 to 60 or 62 to 67
        ranges = []

        c_token = start
        while c_token < len(tokens):
            # note: unwrap with 'degree' converts from degrees to radian and puts it in the -pi to pi range automatically
            lower_bound = tokens.unwrap(c_token, "degree")
            upper_bound = tokens.unwrap(c_token+2, "degree")
            to = tokens.unwrap(c_token+1, "word")
            if c_token+3 < len(tokens):
                or_ = tokens.unwrap(c_token+3, "word")
                if or_ != "or":
                    raise TokenParseException(
                        tokens[c_token+3],
                        "Different allowed angle ranges should be delimited by 'or'."
                    )
            if to != "to":
                raise TokenParseException(
                    tokens[c_token+1],
                    "Angle ranges should be specified with a lower and upper bounds separated by keyword 'to'."
                )
            # check for overlaps
            for other_lower, other_upper in ranges:
                if other_lower < lower_bound < other_upper:
                    raise TokenParseException(
                        tokens[c_token],
                        f"Lower bound falls between an already allowed range {other_lower*180./pi:.2f} to {other_upper*180./pi:.2f}."
                    )
                if other_lower < upper_bound < other_upper:
                    raise TokenParseException(
                        tokens[c_token+2],
                        f"Upper bound falls between an already allowed range {other_lower*180./pi:.2f} to {other_upper*180./pi:.2f}."
                    )

            if not wrap:
                # check for upper < lower
                if upper_bound < lower_bound:
                        raise TokenParseException(
                            tokens[c_token+2],
                            f"Angle upper bound {upper_bound*180./pi:.2f} lower than lower bound {lower_bound*180./pi:.2f}."
                        )
                if lower_bound < 0.:
                    raise TokenParseException(
                        tokens[c_token],
                        f"Angle lower bound must be 0 or larger. Got {lower_bound*180./pi:.2f} instead."
                    )
                ranges.append((lower_bound, upper_bound))

            else:
                # if it's wrapping, pass
                if upper_bound < lower_bound:
                    ranges.append((upper_bound, pi))
                    ranges.append((-pi, lower_bound))
                else:
                    ranges.append((lower_bound, upper_bound))

            if c_token+3 < len(tokens):
                c_token += 4
            else:
                c_token += 3


        # we don't necessarily expect them to be sorted
        ranges.sort()

        # inverting ranges - list of disallowed places
        res = []

        # disallow 0 to first allowed range
        if not wrap and ranges[0][0] > 0.:
            res.append((0., ranges[0][0]))
        elif wrap and ranges[0][0] > -pi:
            res.append((-pi, ranges[0][0]))

        # disallow in between allowed ranges - note, we make sure they don't overlap
        for i in range(len(ranges)-1):
            res.append((ranges[i][1], ranges[i+1][0]))

        if ranges[-1][1] < pi:
            res.append((ranges[-1][1], pi))

        return res

    def line(self, tokens: TokenList) -> None:
        last_reaction = self.parent.d_template
        key = tokens.unwrap(0, "word")
        # sets on which angles and dihedrals were already defined
        angles = {
            (a, b, c, d, e, f)
            for a, b, c, d, e, f, _, _ in last_reaction.angle_limits
        }
        dihedrals = {
            (a, b, c, d, e, f, g, h)
            for a, b, c, d, e, f, g, h, _, _ in last_reaction.dihedral_limits
        }

        possible_keys = {
            "r_max", "r_min", "angle", "dihedral", "rate", "probability"
        }
        if key not in possible_keys:
            close_matches = difflib.get_close_matches(
                key, possible_keys, 1
            )
            raise TokenParseException(
                tokens[0],
                f"Keyword {key} not recognized." +
                (
                    f" Perhaps you meant {close_matches[0]}?"
                    if len(close_matches) > 0 else ""
                ) + f" Valid keywords are: {', '.join(possible_keys)}."
            )

        match key:
            case "r_max":
                idi, atom_i = self.parent.parse_pair(tokens, 1)
                idj, atom_j = self.parent.parse_pair(tokens, 2)
                cutoff = tokens.unwrap(3, "positive")
                last_reaction.add_distance_max((idi, atom_i, idj, atom_j, cutoff))
            case "r_min":
                idi, atom_i = self.parent.parse_pair(tokens, 1)
                idj, atom_j = self.parent.parse_pair(tokens, 2)
                cutoff = tokens.unwrap(3, "positive")
                last_reaction.add_distance_min((idi, atom_i, idj, atom_j, cutoff))
            case "angle":
                idi, atom_i = self.parent.parse_pair(tokens, 1)
                idj, atom_j = self.parent.parse_pair(tokens, 2)
                idk, atom_k = self.parent.parse_pair(tokens, 3)
                if (idi, atom_i, idj, atom_j, idk, atom_k) in angles:
                    raise ParseException(
                        "There is already angle conditions defined for"
                        f" the atoms {(idi, atom_i, idj, atom_j, idk, atom_k)}. "
                        "Please put all angles for these atoms on a single line."
                    )
                for (from_, to) in self.__parse_angle_conditions(tokens, 4, False):
                    last_reaction.add_angle_limit(
                        (idi, atom_i, idj, atom_j, idk, atom_k, cos(from_), cos(to)),
                    )
            case "dihedral":
                idi, atom_i = self.parent.parse_pair(tokens, 1)
                idj, atom_j = self.parent.parse_pair(tokens, 2)
                idk, atom_k = self.parent.parse_pair(tokens, 3)
                idl, atom_l = self.parent.parse_pair(tokens, 4)
                if (idi, atom_i, idj, atom_j, idk, atom_k, idl, atom_l) in dihedrals:
                    raise ParseException(
                        "There is already dihedral conditions defined for"
                        f" the atoms {(idi, atom_i, idj, atom_j, idk, atom_k, idl, atom_l)}. "
                        "Please put all angles for these atoms on a single line."
                    )
                for (from_, to) in self.__parse_angle_conditions(tokens, 5, True):
                    last_reaction.add_dihedral_limit(
                        (idi, atom_i, idj, atom_j, idk, atom_k, idl, atom_l, from_, to),
                    )
            case "rate":
                last_reaction.relative_rate = tokens.unwrap(1, "positive")
            case "probability":
                last_reaction.probability = tokens.unwrap(1, "positive")

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
        return isinstance(parent, ReactionDirective)

    @classmethod
    def get_name(cls) -> str:
        return "conditions"

@register_directive
class BreakDirective(Directive):
    def line(self, tokens: TokenList) -> None:
        self.parent.molecule_type.break_groups.append([
            self.parent.parse_index(tokens, i) for i in range(len(tokens))
        ])

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
        return isinstance(parent, ReactionDirective)

    @classmethod
    def get_name(cls) -> str:
        return "break"


@register_directive
class UpdateDirective(Directive):
    def line(self, tokens: TokenList) -> None:
        self.parent.molecule_type.update_groups.append([
            self.parent.parse_index(tokens, i) for i in range(len(tokens))
        ])

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
        return isinstance(parent, ReactionDirective)

    @classmethod
    def get_name(cls) -> str:
        return "update"


@register_directive
class RedefineDirective(Directive):

    def line(self, tokens: TokenList) -> None:
        reactant_index, atom_index = self.parent.parse_index(tokens, 0, "pair")
        changes = set()
        i = 1
        while i < len(tokens):
            word = tokens.unwrap(i, "word")
            if word in changes:
                raise TokenParseException(
                    tokens[i],
                    f"Duplicate entry {word}."
                )
            changes.add(word)
            match word:
                case "name":
                    self.parent.molecule_type.renames.append((
                        reactant_index, atom_index, tokens.unwrap(i+1, "word")
                    ))
                case "type":
                    self.parent.molecule_type.retypes.append((
                        reactant_index, atom_index, tokens.unwrap(i+1, "word")
                    ))
                case "charge":
                    self.parent.molecule_type.recharges.append((
                        reactant_index, atom_index, tokens.unwrap(i+1, "float")
                    ))
                case "mass":
                    self.parent.molecule_type.remasses.append((
                        reactant_index, atom_index, tokens.unwrap(i+1, "float")
                    ))
                case _:
                    raise TokenParseException(
                        tokens[1 + i * 2],
                        f"Unknown atom property {word}."
                    )

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
        return isinstance(parent, ReactionDirective)

    @classmethod
    def get_name(cls) -> str:
        return "redefine"

@register_directive
class SoftCoreDirective(Directive):

    def line(self, tokens: TokenList) -> None:
        reactant_index, atom_index = self.parent.parse_index(tokens, 0, "pair")
        sc_lam = tokens.unwrap(1, "float")
        sc_alpha = tokens.unwrap(2, "float")
        self.parent.molecule_type.soft_core.append((reactant_index, atom_index, sc_lam, sc_alpha))

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
        return isinstance(parent, ReactionDirective)

    @classmethod
    def get_name(cls) -> str:
        return "softcore"
