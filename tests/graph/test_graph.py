import os
import pytest
import json
from martini_daemon.topstar import TopStar
from martini_daemon.top_parser import DaemonTopFile
from martini_daemon.forces.nonbonded import NonBonded
from martini_daemon.fragment import Fragment


# == CONFIG ==
@pytest.fixture
def rootdir(request):
    return os.path.dirname(request.path)


tests = [
    "single", "square", "star", "BDT", "opt", "equiv", "v_shape", "bicycle",
    "spiked_triangle"
]


# == TEST CLASS ==
class TestGraph():
    def get_topology(self, path: str) -> TopStar:
        assert os.path.isfile(path)
        ok, res = DaemonTopFile(path, NonBonded())
        assert ok
        _, top = res
        return top

    def parse_expected(self, path) -> list[(str, list[int])]:
        """
            Parses a json of the format
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
            data = [
                (line[0], [x - 1 for x in line[1:]]) for line in json.load(f)
            ]
        return data

    def try_match(self, frag: Fragment, name: str, parts: list[int]):
        if frag.name != name:
            return False
        frag_parts = list(filter(lambda x: x != -1, frag.atoms))
        if len(frag_parts) != len(parts):
            return False
        for i, atom in enumerate(frag_parts):
            if atom != parts[i]:
                return False
        return True

    def compare(self, topstar: TopStar, expected: list[(str, list[int])]):
        """
        Makes sure topstar has all of and only the fragments in expected.
        Expected contains the frag names and particles.
        Note: the order of expected is arbitrary.
        """
        matched_frags = set()
        for (frag_name, part_ids) in expected:
            frag_ids = topstar.defrag_list[part_ids[0]]
            frags = [topstar.frag_list[frag_id] for frag_id in frag_ids]
            matched = False
            for frag in frags:
                if frag in matched_frags:
                    continue
                if self.try_match(frag, frag_name, part_ids):
                    matched = True
                    matched_frags.add(frag)
                    break
            if not matched:
                for frag in frags:
                    print("found frags:")
                    # convert to 1 based indexing
                    atoms = [
                        atom+1 if atom >= 0 else atom for atom in frag.atoms
                    ]
                    print(f"frag {frag.name} {atoms}")
                assert False
        assert len(expected) == len(topstar.frag_list)

    @pytest.mark.parametrize("x", tests)
    def test_graph(self, x, rootdir):
        os.chdir(rootdir)
        assert os.path.isdir(x)
        os.chdir(x)
        topstar = self.get_topology("system.top")
        expected = self.parse_expected("result.json")
        self.compare(topstar, expected)
