from ..utils import backup_try
import zlib


class Reporter():
    """
    Placeholder / base class for DaemonSimulation reporters.
    """

    # filled out by DaemonSimulation
    _topstar = None
    _sysstar = None

    # ====== (Optional) Helpers ======
    _handle = None
    _comp = None
    _mode = None

    def _open(self, path):
        """
        Set self._handle for an uncompressed text format.
        """
        assert self._handle is None and self._comp is None, "_open must only be called once."
        backup_try(path)
        self._handle = open(path, "w")
        self._mode = "w"
        self._comp = None

    def _open_compressed(self, path):
        """
        Set self._handle for a compressed binary format
        (with self._comp for compression).
        """
        assert self._handle is None and self._comp is None, "_open must only be called once."
        backup_try(path)
        self._handle = open(path, "wb")
        self._mode = "wb"
        self._comp = zlib.compressobj(6)

    def _write(self, a):
        if self._comp is not None:
            self._handle.write(self._comp.compress(a))
        else:
            self._handle.write(a)
        self._handle.flush()

    def _print(self, a):
        assert self._comp is None and self._mode == "w", "only use _print with uncompressed text modes"
        self._handle.write(a + "\n")
        self._handle.flush()

    def __del__(self):
        self.finish()

    def finish(self):
        if self._handle is not None:
            if self._comp is not None:
                self._handle.write(self._comp.flush())
                self._comp = None
            self._handle.close()
            self._handle = None

    # ====== S* reporting overloadable ======
    # on_set_xtc_path, on_xtc_frame, on_write_gro -> should be used to reporters that try to mirror trajectories
    def on_set_xtc_path(self, xtc_name):
        """
            Called when xtc_name is first known, before writing any xtc.

            Useful to run backups on any file a reporter writes that tries
            to parallel trajectories, as set_xtc_path in S* backs up .xtc
            files so it can start a fresh one.

            Do not rely on this to initialize! May not be called when loading
            from checkpoints.
        """
        pass

    def on_xtc_frame(self, i, pos, box, xtc_name):
        """
            Called every time after the XTC is written to
            pos is the positions passed to the XTC, name is the name of the xtc
            without the .xtc extension
        """
        pass

    def on_write_gro(self, pos, box, gro_name):
        # Called after any .gro coordinates file is written

        # pos is the positions written to the GRO file, gro_name is the name
        # of the gro file without the .gro extension
        pass

    # ====== T* reporting overloadable ======
    # init_dm, pre_detection, pre_modification, post_modification -> should be used for reporters that want to report on T* and D/M
    def init_dm(self, name):
        pass

    def pre_detection(self, i, name):
        # called before the D algorithm runs
        pass

    def pre_modification(self, i, reactions, name):
        # called before the modification algorithm with a list of reactions
        # to be processed in the modification algorithm
        # only called if len(reactions) > 0
        # beware: some reactions can be rejected by the modification algorithm
        # use post_modification for accurate reactions
        pass

    def post_modification(self, i, name, reactions):
        # called after the modification algorithm, before the system is
        # reinitialized. Only runs if there were reactions.
        pass

    # ====== DaemonIntegrator hooks ======

    def post_reinitialize(self, i):
        # called after reinitializing post-reactions
        pass

    def post_sc_enable(self, i):
        # called after enabling soft core
        pass

    def post_di_minimize(self):
        # called after daemon integrator minimized
        pass

    def post_di_equilibrate(self):
        # called after daemon integrator has equilibrated
        pass

    def post_sc_disable(self, i):
        # called after soft core is disabled
        pass


    # ====== Misc overloadable ======
    # adds something to the line shown interactively to the user
    def interactive_line(self) -> str:
        pass


def read_compressed(path):
    with open(path, "rb") as f:
        return zlib.decompress(f.read())
