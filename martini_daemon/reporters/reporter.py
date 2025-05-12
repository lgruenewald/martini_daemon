class Reporter():
    """Placeholder / base class for reporters. _sysstar and _topstar are
    added as fields when a reporter is added to a DaemonSimulation."""

    # ====== S* reporting ======
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

    # ====== T* reporting ======
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
        pass

    def post_modification(self, i, name):
        # called after the modification algorithm, before the system is
        # reinitialized. Only runs if there were reactions.
        pass

    # ====== Misc ======
    # adds something to the line shown interactively to the user
    def interactive_line(self) -> str:
        pass
