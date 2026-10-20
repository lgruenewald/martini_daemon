"""Test the internal function _parse_angle_conditions."""

import math

from martini_daemon import TokenList
from martini_daemon.__topstar.reaction_directives import _parse_angle_conditions


def test_parse_angle_conditions() -> None:
    """Test the internal function _parse_angle_conditions."""
    lines = [
        (False, "1.0 to 50", [(0.0, 1.0), (50, 180)]),
        (False, "8.0 to 140 or 150 to 156", [(0, 8), (140, 150), (156, 180)]),
        (True, "8.0 to 140 or 150 to 156", [(-180, 8), (140, 150), (156, 180)]),
        (
            True,
            "-150 to -130 or 0 to 20 or 170 to -175",
            [(-175, -150), (-130, 0), (20, 170)],
        ),
        (True, "170 to -170", [(-170, 170)]),
    ]
    for wrap, line, exp in lines:
        tl = TokenList(line, "test", 1)
        res = _parse_angle_conditions(tl, 0, wrap)
        for (res_s, res_e), (exp_s, exp_e) in zip(res, exp):
            exp_s *= math.pi / 180.0
            exp_e *= math.pi / 180.0
            assert math.isclose(res_s, exp_s), f"res_s {res_s}, exp_s {exp_s}"
            assert math.isclose(res_e, exp_e), f"res_e {res_e}, exp_e {exp_e}"
