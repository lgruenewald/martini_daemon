from .reporter import Reporter
from ..utils import backup_try, cross_box
import zlib


class VMDBondReporter(Reporter):
    """Reporter that prints a live list of bonds to a file, that can be
    visualized in VMD using daemon.tcl

    Limitation right now:
    Does not report bonds that cross the periodic boundary conditions,
    once a daemon aware pbc whole is written perhaps this can change"""

    handle = None

    def __init__(self):
        pass

    def on_set_xtc_path(self, xtc_name):
        self.bonds_name = xtc_name + ".bonds"
        backup_try(self.bonds_name)
        self.comp_obj = zlib.compressobj(6)
        self.handle = open(self.bonds_name, "wb")

    def on_xtc_frame(self, i, pos, box, _):
        bonds = "{" + "} {".join(
            map(
                lambda part:  # part[0] - i, part[1] - list of interactions
                " ".join(
                    map(
                        lambda ij:  # ij[0] - i, ij[1] - j
                            str(ij[1]) if ij[0] == part[0] else str(ij[0]),
                        map(
                            lambda inter: inter.get_members(),
                            filter(
                                lambda inter:
                                    not cross_box(
                                        pos[inter.get_members()[0]],
                                        pos[inter.get_members()[1]],
                                        box
                                    ),
                                filter(
                                    lambda inter: inter.is_instance("bond"),
                                    part[1]
                                )
                            )
                        )
                    )
                ),
                enumerate(self._topstar.interaction_list)
            )
        ) + "} "

        flat = self.comp_obj.compress(bonds.encode("utf-8"))

        self.handle.write(flat)

    def __del__(self):
        if self.handle is not None:
            self.handle.write(self.comp_obj.flush())
            self.handle.close()
