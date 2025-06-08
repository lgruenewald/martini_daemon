import os
import glob
import pytest
from martini_daemon.simulation import Simulation
from martini_daemon.reporters.topstar import ReactionReporter


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
    "distance",
    "angle",
# TODO
#    "dihedral",
#    "mono_distance",
#    "optional",
#    "overlap"
]


# == TEST CLASS ==
class TestDetection():

    def get_sim(self, top, gro):
        rep = ReactionReporter()
        sim = Simulation(
            top, gro,
            reporters=[
                rep
            ]
        )
        sim.step(0, xtc=False, dm=True)
        # force closing of file
        # TODO oof
        rep.__del__()
        return "out.reactions"

    def read_reactions(self, path):
        """
        .reactions format reader suited for test_detection.py
        """
        reactions = []
        with open(path, "r") as f:
            lines = f.read().splitlines()
            for line in lines:
                if line[0] == "#" or len(line) == 0:
                    continue
                elems = line.split(";")
                frame, rx = elems[0].split(",")
                # list of atoms
                frags = [
                    elem.split(",")[2:] for elem in elems[1:]
                ]
                reactions.append(
                    (rx, frags)
                )
        return reactions

    def compare(self, reactions, expected):
        dump = f"\nGot: {reactions}, expected: {expected}."
        assert len(reactions) == len(expected), f"first len check {dump}"
        for (r1, frags1), (r2, frags2) in zip(reactions, expected):
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

    @pytest.mark.parametrize("x", tests)
    def test_detection(self, x, rootdir):
        # enter dir
        os.chdir(rootdir)
        assert os.path.isdir(x)
        os.chdir(x)

        # meat of the test
        rx_path = self.get_sim("system.top", "system.gro")
        reactions = self.read_reactions(rx_path)
        expected = self.read_reactions("expected.reactions")
        self.compare(reactions, expected)

        # cleanup
        for filename in glob.glob("./out*"):
            os.remove(filename)
        for filename in glob.glob("./#*"):
            os.remove(filename)
