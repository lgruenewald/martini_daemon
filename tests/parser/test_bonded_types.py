"""Test [bondtypes], [angletypes] and [dihedraltypes] parameter lookup."""

import math
import os
from pathlib import Path

import pytest

import martini_daemon  # noqa: F401 (registers directives and forces)
from martini_daemon.__core import System
from martini_daemon.__parser.gromacs_top_file import (
    GromacsTopFile,
    InvalidTopologyError,
)

ITP = os.path.join(os.path.dirname(__file__), "..", "..", "itp", "martini_v3.0.0.itp")

# four atoms with distinct types, so a lookup on the wrong atoms shows up
ATOMS = """
[ atoms ]
1 P1 1 MOL A 1 0.0 72.0
2 P2 1 MOL B 2 0.0 72.0
3 P3 1 MOL C 3 0.0 72.0
4 P4 1 MOL D 4 0.0 72.0
"""


def parse(
    tmp_path: Path, types: str, interactions: str, empty_types: bool = True
) -> list[tuple[str, list[int], list[float], bool]]:
    """Parse a one-molecule topology and return the molecule's interactions."""
    if empty_types:
        # empty headers so the lookup tests don't depend on whether
        # the *types directives are mandatory
        types = "[ bondtypes ]\n[ angletypes ]\n[ dihedraltypes ]\n" + types
    top = tmp_path / "system.top"
    top.write_text(
        f'#include "{os.path.abspath(ITP)}"\n'
        f"{types}\n"
        "[ moleculetype ]\nMOL 1\n"
        f"{ATOMS}\n{interactions}\n"
        "[ system ]\ntest\n\n[ molecules ]\nMOL 1\n"
    )
    system = System()
    GromacsTopFile(system, str(top))
    return system.molecule_types["MOL"].interactions


def params(
    interactions: list[tuple[str, list[int], list[int | float], bool]],
) -> list[list[float]]:
    """Extract params from moleculetype object .interactions field."""
    return [p for _, _, p, _ in interactions]


# === basics ===


def test_no_types_directives_needed(tmp_path: Path) -> None:
    """Topologies without any *types directive must still parse."""
    got = parse(tmp_path, "", "[ bonds ]\n1 2 1 0.3 1000\n", empty_types=False)
    assert params(got) == [[0.3, 1000.0]]


def test_bond_lookup(tmp_path: Path) -> None:
    got = parse(
        tmp_path,
        "[ bondtypes ]\nP1 P2 1 0.47 1250\n",
        "[ bonds ]\n1 2 1\n",
    )
    assert params(got) == [[0.47, 1250.0]]


def test_bond_lookup_uses_the_right_atoms(tmp_path: Path) -> None:
    """Bond 2-3 must use the P2 P3 entry, not the P1 P2 one."""
    got = parse(
        tmp_path,
        "[ bondtypes ]\nP1 P2 1 0.1 100\nP2 P3 1 0.2 200\nP3 P4 1 0.3 300\n",
        "[ bonds ]\n2 3 1\n",
    )
    assert params(got) == [[0.2, 200.0]]


def test_bond_lookup_reversed(tmp_path: Path) -> None:
    got = parse(
        tmp_path,
        "[ bondtypes ]\nP2 P1 1 0.47 1250\n",
        "[ bonds ]\n1 2 1\n",
    )
    assert params(got) == [[0.47, 1250.0]]


def test_explicit_params_override(tmp_path: Path) -> None:
    got = parse(
        tmp_path,
        "[ bondtypes ]\nP1 P2 1 0.47 1250\n",
        "[ bonds ]\n1 2 1 0.3 1000\n",
    )
    assert params(got) == [[0.3, 1000.0]]


def test_lookup_is_per_line(tmp_path: Path) -> None:
    """Each bond gets its own parameters."""
    got = parse(
        tmp_path,
        "[ bondtypes ]\nP1 P2 1 0.1 100\nP2 P3 1 0.2 200\n",
        "[ bonds ]\n1 2 1\n2 3 1\n1 2 1\n",
    )
    assert params(got) == [[0.1, 100.0], [0.2, 200.0], [0.1, 100.0]]


def test_parameterless_type_needs_no_lookup(tmp_path: Path) -> None:
    """Bond type 5 has no parameters, so 1 2 5 is complete as is."""
    got = parse(tmp_path, "", "[ bonds ]\n1 2 5\n")
    assert params(got) == [[]]


def test_missing_type_fails(tmp_path: Path) -> None:
    with pytest.raises(InvalidTopologyError):
        parse(tmp_path, "[ bondtypes ]\nP1 P2 1 0.47 1250\n", "[ bonds ]\n3 4 1\n")


def test_function_type_must_match(tmp_path: Path) -> None:
    """A func 1 entry must not be used for a func 2 bond."""
    with pytest.raises(InvalidTopologyError):
        parse(tmp_path, "[ bondtypes ]\nP1 P2 1 0.47 1250\n", "[ bonds ]\n1 2 2\n")


@pytest.mark.parametrize("second", ["P1 P2", "P2 P1"])
def test_duplicate_fails(tmp_path, second) -> None:
    """Defining the same entry twice (in either order) is an error."""
    with pytest.raises(InvalidTopologyError):
        parse(
            tmp_path,
            f"[ bondtypes ]\nP1 P2 1 0.1 100\n{second} 1 0.2 200\n",
            "[ bonds ]\n1 2 1\n",
        )


def test_duplicate_other_func_ok(tmp_path: Path) -> None:
    """The same atom types with a different function type are not duplicates."""
    got = parse(
        tmp_path,
        "[ bondtypes ]\nP1 P2 1 0.1 100\nP1 P2 2 0.2 200\n",
        "[ bonds ]\n1 2 2\n",
    )
    assert params(got) == [[0.2, 200.0]]


def test_symmetric_entry_is_not_a_duplicate(tmp_path: Path) -> None:
    """P1 P1 reversed is itself, which must not count as a duplicate."""
    got = parse(
        tmp_path,
        "[ bondtypes ]\nP1 P1 1 0.1 100\nX P2 1 0.2 200\n",
        "[ bonds ]\n1 2 1\n",
    )
    assert params(got) == [[0.2, 200.0]]


def test_define_in_type_params(tmp_path: Path) -> None:
    got = parse(
        tmp_path,
        "#define B_P1P2 0.47 1250\n[ bondtypes ]\nP1 P2 1 B_P1P2\n",
        "[ bonds ]\n1 2 1\n",
    )
    assert params(got) == [[0.47, 1250.0]]


# === angles ===


def test_angle_lookup(tmp_path: Path) -> None:
    got = parse(
        tmp_path,
        "[ angletypes ]\nP1 P2 P3 2 120 25\n",
        "[ angles ]\n1 2 3 2\n",
    )
    assert params(got) == [[pytest.approx(math.radians(120)), 25.0]]


def test_angle_lookup_reversed(tmp_path: Path) -> None:
    got = parse(
        tmp_path,
        "[ angletypes ]\nP3 P2 P1 2 120 25\n",
        "[ angles ]\n1 2 3 2\n",
    )
    assert params(got) == [[pytest.approx(math.radians(120)), 25.0]]


def test_angle_only_full_reversal(tmp_path: Path) -> None:
    """P2 P1 P3 is a different angle from P1 P2 P3."""
    with pytest.raises(InvalidTopologyError):
        parse(
            tmp_path,
            "[ angletypes ]\nP2 P1 P3 2 120 25\n",
            "[ angles ]\n1 2 3 2\n",
        )


# === dihedrals ===


def test_dihedral_lookup(tmp_path: Path) -> None:
    got = parse(
        tmp_path,
        "[ dihedraltypes ]\nP1 P2 P3 P4 1 180 10 2\n",
        "[ dihedrals ]\n1 2 3 4 1\n",
    )
    assert params(got) == [[pytest.approx(math.pi), 10.0, 2.0]]


def test_dihedral_lookup_reversed(tmp_path: Path) -> None:
    got = parse(
        tmp_path,
        "[ dihedraltypes ]\nP4 P3 P2 P1 1 180 10 2\n",
        "[ dihedrals ]\n1 2 3 4 1\n",
    )
    assert params(got) == [[pytest.approx(math.pi), 10.0, 2.0]]


def test_dihedral_outer_wildcards(tmp_path: Path) -> None:
    got = parse(
        tmp_path,
        "[ dihedraltypes ]\nX P2 P3 X 1 180 10 2\n",
        "[ dihedrals ]\n1 2 3 4 1\n",
    )
    assert params(got) == [[pytest.approx(math.pi), 10.0, 2.0]]


def test_dihedral_inner_wildcards(tmp_path: Path) -> None:
    """AMBER-style improper with wildcards in positions 1 and 2."""
    got = parse(
        tmp_path,
        "[ dihedraltypes ]\nX X P3 P4 1 180 10 2\n",
        "[ dihedrals ]\n1 2 3 4 1\n",
    )
    assert params(got) == [[pytest.approx(math.pi), 10.0, 2.0]]


@pytest.mark.parametrize("specific_first", [True, False])
def test_dihedral_most_specific_wins(tmp_path: Path, specific_first: bool) -> None:
    """Fewest wildcards wins, regardless of order in the file."""
    specific = "P1 P2 P3 P4 1 0 1 1\n"
    generic = "X P2 P3 X 1 180 10 2\n"
    lines = specific + generic if specific_first else generic + specific
    got = parse(tmp_path, "[ dihedraltypes ]\n" + lines, "[ dihedrals ]\n1 2 3 4 1\n")
    assert params(got) == [[0.0, 1.0, 1.0]]


@pytest.mark.parametrize("one_first", [True, False])
def test_dihedral_fewer_wildcards_wins(tmp_path: Path, one_first: bool) -> None:
    one = "P1 P2 P3 X 1 0 1 1\n"
    two = "X P2 P3 X 1 180 10 2\n"
    lines = one + two if one_first else two + one
    got = parse(tmp_path, "[ dihedraltypes ]\n" + lines, "[ dihedrals ]\n1 2 3 4 1\n")
    assert params(got) == [[0.0, 1.0, 1.0]]


def test_dihedral_wildcard_reversed(tmp_path: Path) -> None:
    """A wildcard entry written in reverse order still matches."""
    got = parse(
        tmp_path,
        "[ dihedraltypes ]\nP4 P3 X X 1 180 10 2\n",
        "[ dihedrals ]\n1 2 3 4 1\n",
    )
    assert params(got) == [[pytest.approx(math.pi), 10.0, 2.0]]


# === wildcards in bonds/angles (allowed by choice) ===


def test_bond_wildcard(tmp_path: Path) -> None:
    got = parse(
        tmp_path,
        "[ bondtypes ]\nX P2 1 0.47 1250\n",
        "[ bonds ]\n1 2 1\n",
    )
    assert params(got) == [[0.47, 1250.0]]
