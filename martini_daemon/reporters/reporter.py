class Reporter():

    def __init__(self, sysstar, topstar):
        self._sysstar = sysstar
        self._topstar = topstar

    def on_set_xtc_path(self, xtc_name):
        # called when xtc_name is first known, before writing any xtc
        pass

    def on_xtc_frame(self, pos, box, xtc_name):
        # Called every time after the XTC is written to
        # pos is the positions passed to the XTC, name is the name of the xtc
        # without the .xtc extension
        pass

    def on_write_gro(self, pos, box, gro_name):
        # Called after a GRO is written
        # pos is the positions written to the GRO file, gro_name is the name
        # of the gro file without the .gro extension
        pass

    def pre_detection(self):
        # called before the D algorithm runs
        pass

    def pre_modification(self, reactions):
        # called before the modification algorithm with a list of reactions
        # to be processed in the modification algorithm
        # only called if len(reactions) > 0
        pass

    def post_modification(self):
        # called after the modification algorithm, before the system is
        # reinitialized. Only runs if there were reactions.
        pass
