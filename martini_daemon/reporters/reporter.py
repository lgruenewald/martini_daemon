class Reporter():

    def __init__(self, sysstar):
        self._sysstar = sysstar

    def pre_steps(self, xtc_name):
        # called when xtc_name is first known, before writing any xtc
        pass

    def step(self, pos, box, xtc_name):
        # Called every time after the XTC is written to
        # pos is the positions passed to the XTC, name is the name of the xtc
        # without the .xtc extension
        pass

    def final(self, pos, box, gro_name):
        # Called after a GRO is written
        # pos is the positions written to the GRO file, gro_name is the name
        # of the gro file without the .gro extension
        pass
