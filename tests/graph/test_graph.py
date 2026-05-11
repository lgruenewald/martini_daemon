"""Test the graph matching algorithm."""

import json
import os

import pytest

from martini_daemon import Fragment, Simulation, TopStar


# == CONFIG ==
@pytest.fixture
def rootdir(request: pytest.FixtureRequest) -> str:
    """Return the root directory for this test."""
    return os.path.dirname(request.path)


tests = [
    "single",
    "square",
    "star",
    "BDT",
    "opt",
    "equiv",
    "v_shape",
    "bicycle",
    "spiked_triangle",
]


# == TEST CLASS ==
def get_topology(path: str) -> TopStar:
    """Return a TopStar instance from a topology file."""
    sim = Simulation(path, None, 0, [])
    sim.finish()
    return sim.top


def parse_expected(path: str) -> list[tuple[str, list[int]]]:
    """
    Parse a json of excepted graph match results.

    Format:
    [
        ["name1", *parts],
        ...
        ["namen", *parts]
    ]
    where *parts is a list of atom indices (1 indexed).
    Converts the indexing to 0 based.
    """
    assert os.path.isfile(path)
    with open(path) as f:
        return [(line[0], [x - 1 for x in line[1:]]) for line in json.load(f)]


def try_match(frag: Fragment, name: str, parts: list[int]) -> bool:
    """Check if a fragment matches the expected graph."""
    if frag.name != name:
        return False
    frag_parts = list(filter(lambda x: x != -1, frag.atoms))
    if len(frag_parts) != len(parts):
        return False
    return all(atom == parts[i] for i, atom in enumerate(frag_parts))


def print_error(
    name: str, frags: list[Fragment | None], expected: list[tuple[str, list[int]]]
) -> None:
    """Print error message and details of frags and expected."""
    print(f"found frags for {name}:")
    # convert to 1 based indexing
    for frag in frags:
        assert frag is not None
        atoms = [atom + 1 if atom >= 0 else atom for atom in frag.atoms]
        print(f"frag {frag.name} {atoms}")
    print("expected frags:")
    for exp in expected:
        print(exp)


def compare(name: str, top: TopStar, expected: list[tuple[str, list[int]]]) -> None:
    """
    Make sure top has all of and only the fragments in expected.

    Expected contains the frag names and particles.
    Note: the order of expected is arbitrary.
    """
    matched_frags = set()
    all_frags = []
    for frag_name, part_ids in expected:
        frag_ids = top.frag_list.frag_ids_for(part_ids[0])
        frags = [top.frag_list.get_fragment(frag_id) for frag_id in frag_ids]
        all_frags += frags
        matched = False
        for frag in frags:
            if frag in matched_frags:
                continue
            if try_match(frag, frag_name, part_ids):
                matched = True
                matched_frags.add(frag)
                break
        if not matched:
            print_error(name, frags, expected)
            assert False
    if len(expected) != top.frag_list.num_fragments():
        print_error(name, all_frags, expected)
        assert False


@pytest.mark.parametrize("x", tests)
def test_graph(x: str, rootdir: str) -> None:
    """Test the graph matching algorithm by comparing results with expected .json files."""
    os.chdir(rootdir)
    assert os.path.isdir(x)
    os.chdir(x)
    top = get_topology("system.top")
    expected = parse_expected("result.json")
    compare(x, top, expected)
