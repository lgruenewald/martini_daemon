from ..old_sysstar import SysStar
from ..old_topstar import TopStar
from .reporter import Reporter
import pickle
import os

# TODO needs automated testing
# TODO save() and load() should be part of constructors
# TODO should this be removed?
class CheckpointReporter(Reporter):
    """
        Reporter that enables simple continuable
        simulations, by writing a checkpoint (.chk)
        file every once in a while.

        Checkpoints can be used to restore DaemonSimulation
        instances with the same T*, S*, openmm system
        as was saved during checkpoint writing.

        the openmm part of checkpoint writing is
        inspired by the their checkpoint reporter
        described at https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/checkpointreporter.py

        every_n_frames gives a goal checkpoint frequency (in MD steps).
        The frequency might not be perfect -- it is only checked
        on xtc frame writes, so every_n_frames gives
        a lower bound on the frequency of writing
        checkpoints.
    """

    def __init__(self, every_n_frames):
        self.every_n_frames = every_n_frames
        self.last_frame = 0
        self.since_checkpoint = 0

    def on_xtc_frame(self, i, pos, box, xtc_name):
        self.since_checkpoint += i - self.last_frame
        self.last_frame = i
        if i == 0 or self.since_checkpoint < self.every_n_frames:
            return
        check_path = xtc_name + ".chk"
        save_checkpoint(check_path, i, box, self._topstar, self._sysstar)
        self.since_checkpoint = 0


def load_checkpoint(
    path, logger, nonbonded_force, nlist_cutoff, max_abs_rate, highest_prob
):
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
        sys = SysStar(logger, nonbonded_force)
        top = TopStar(
            sys, logger, nlist_cutoff, max_abs_rate, highest_prob, None
        )
        sys.load(f)
        top.load(f)
        assert f.load() == "magic"
        return i, sys, top


def save_checkpoint(path, i, box, top: TopStar, sys: SysStar):
    """
        Atomically saves a checkpoint file at <path>.

        Saves the current frame (i), T* (top), S* (sys).

        Load with load_checkpoint.
    """
    # power outage or similar during checkpoint writing?
    # we need atomic writes, we achieve this with
    # writing a tempfile first, and then using os.replace
    # to move the file, overwriting the old one
    tmpfile = f"{path}#"
    with open(tmpfile, "wb") as file:
        f = pickle.Pickler(file)
        f.dump(i)
        sys.save(f, box)
        top.save(f)
        f.dump("magic")
    os.replace(tmpfile, path)
