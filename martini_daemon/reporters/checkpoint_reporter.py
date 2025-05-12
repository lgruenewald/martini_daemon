# the openmm part of reporting is inspired by
# https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/checkpointreporter.py

from ..sysstar import SysStar
from ..topstar import TopStar
from .reporter import Reporter
import pickle


class CheckpointReporter(Reporter):

    def __init__(self):
        pass

    def on_xtc_frame(self, i, pos, box, xtc_name):
        if i == 0:
            return
        check_path = xtc_name + ".chk"
        save_checkpoint(check_path, i, box, self._topstar, self._sysstar)


def load_checkpoint(path, logger, epsilon_r, nonbonded_cutoff, nlist_cutoff):
    """
        Loads checkpoint file at <path>.

        Returns:
        - the current frame / progression of the simulation (i)
        - T*
        - S*, including an openmm context with positions and velocities

        This essentially replaces loading a .gro and a .top file,
        but the remaining parameters will stay as defined in
        DaemonSimulation().
    """
    with open(path, "rb") as file:
        f = pickle.Unpickler(file)
        i = f.load()
        sys = SysStar(logger, epsilon_r, nonbonded_cutoff)
        top = TopStar(sys, logger, nlist_cutoff)
        sys.load(f)
        top.load(f)
        assert f.load() == "magic"
        return i, sys, top


def save_checkpoint(path, i, box, top: TopStar, sys: SysStar):
    """
        Saves a checkpoint file at <path>.

        Saves the current frame (i), T* (top), S* (sys).

        Load with load_checkpoint.
    """
    with open(path, "wb") as file:
        f = pickle.Pickler(file)
        f.dump(i)
        sys.save(f, box)
        top.save(f)
        f.dump("magic")
