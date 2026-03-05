import os
import glob
import shutil
import pytest
from martini_daemon.simulation import Simulation
from martini_daemon.reporters.sysstar_dump import SysStarDump
from martini_daemon.reporters.topstar import FragCountReporter, ReactionReporter
from math import isclose


# == CONFIG ==
@pytest.fixture
def rootdir(request):
    return os.path.dirname(request.path)


tests = [
    # adding new interactions
    "bond_make",
    # [break]
    "bond_break", "dihedral_break",
    # [update]
    "update_group", "update_group4"
]


# == TEST CLASS ==
class TestDetection:
    @staticmethod
    def get_sim(top: str, gro: str) -> str:
        rep = SysStarDump()
        sim = Simulation(
            top, gro,
            reporters=[
                rep, FragCountReporter(), ReactionReporter()
            ]
        )
        sim.step(0, xtc=False, dm=True)
        # force closing of file
        # TODO oof
        # it is what it is
        rep.__del__()
        return "out.sstar"

    @staticmethod
    def compare(dump_new: str, dump_reference: str) -> None:
        reference = SysStarDump.read_dump(dump_reference)
        new = SysStarDump.read_dump(dump_new)
        assert len(new) == len(reference)
        # structural equality with floats
        # TODO surely this is an opportunity to discover something cool and new and simple
        # TODO make it order insensitive for both forces and within forces
        for i in range(len(new)):
            frame, atoms, forces = new[i]
            frame_ref, atoms_ref, forces_ref = reference[i]
            assert frame == frame_ref
            for j in range(len(atoms)):
                assert atoms[j][0] == atoms_ref[j][0] # name
                assert atoms[j][1] == atoms_ref[j][1] # resid
                assert atoms[j][2] == atoms_ref[j][2] # res_name
                assert atoms[j][3] == atoms_ref[j][3] # atom_type
                assert isclose(atoms[j][4], atoms_ref[j][4]) # charge
                assert isclose(atoms[j][5], atoms_ref[j][5]) # mass
            assert len(forces) == len(forces_ref)
            # each force obj
            for (name, force), (name_ref, force_ref) in zip(forces, forces_ref):
                assert name == name_ref
                # single line of force
                for inter, inter_ref in zip(force, force_ref):
                    # atom ids and params
                    for val, val_ref in zip(inter, inter_ref):
                        assert isclose(val, val_ref)


    @pytest.mark.parametrize("x", tests)
    def test_detection(self, x, rootdir):
        # enter dir
        os.chdir(rootdir)
        assert os.path.isdir(x)
        os.chdir(x)

        # meat of the test
        dump_new = TestDetection.get_sim("system.top", "system.gro")
        dump_reference = "reference.sstar"
        TestDetection.compare(dump_new, dump_reference)

        # cleanup
        for filename in glob.glob("./out*"):
            os.remove(filename)
        for filename in glob.glob("./#*"):
            os.remove(filename)
