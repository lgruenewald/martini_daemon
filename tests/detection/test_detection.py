import os
import glob
import pytest
from martini_daemon import Simulation, ReactionReporter

# == CONFIG ==
@pytest.fixture
def rootdir(request):
    return os.path.dirname(request.path)


tests = [
    # real systems
    "tetrahedral",
    "polyurethane",
    "bdt",
    # dummy systems
    "three_reagent",
    "distance",
    "angle",
    "dihedral",
    "optional",
    "overlap"
]


def get_sim(top: str, gro: str) -> list[tuple[str, str, list[list[int]]]]:
    rep = ReactionReporter()
    sim = Simulation(
        top, gro, 0,
        reporters=[
            rep
        ]
    )
    sim.step(0, traj=False, dm=True)
    sim.finish()
    return ReactionReporter.read_reactions("out.reactions")

def compare(reactions, expected):
    dump = f"\nGot: {reactions}, expected: {expected}."
    for (_, r1, frags1), (_, r2, frags2) in zip(reactions, expected):
        # while order in theory can be different, it is simpler
        # to for now make systems where we just form the expected
        # .reactions in the order detection will (deterministically)
        # output them
        assert r1 == r2, f"Name differs. {dump}"
        assert len(frags1) == len(frags2), "second len check" + dump
        for atoms1, atoms2 in zip(frags1, frags2):
            assert len(atoms1) == len(atoms2), "third len check" + dump
            for atom1, atom2 in zip(atoms1, atoms2):
                assert atom1 == atom2, f"Atom differs. {dump}"
    assert len(reactions) == len(expected), f"first len check {dump}"

@pytest.mark.parametrize("x", tests)
def test_detection(x, rootdir):
    # enter dir
    os.chdir(rootdir)
    assert os.path.isdir(x)
    os.chdir(x)

    # vegan meat substitute of the test
    reactions = get_sim("system.top", "system.gro")
    expected = ReactionReporter.read_reactions("expected.reactions")
    compare(reactions, expected)

    # cleanup
    for filename in glob.glob("./out*"):
        os.remove(filename)
    for filename in glob.glob("./#*"):
        os.remove(filename)
