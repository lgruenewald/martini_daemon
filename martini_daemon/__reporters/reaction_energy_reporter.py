import os

from ..__reporter import Reporter
from ..__simulation import Simulation
from .variables_reporter import write_energies


class ReactionEnergyReporter(Reporter):
    def __init__(self, write_coords=False, ext=".gro"):
        """Reporter that will write energies before and after a reaction.

        :param write_coords: if set to True, it will print .gro files pre and post minimization
        :param ext: extension for writing the geometries. Set to ".xyz" if xyz files are desired.
        """
        self.write_coords = write_coords
        self.ext = ext

    def __truncate(self, sim: Simulation):
        h = sim.get_handle(".rxener")
        h.seek(0, os.SEEK_SET)
        prev_pos = 0
        while (line := h.readline()) != b"":
            frame = int(line.decode("utf-8").split(",")[1])
            if frame > sim.current_step:
                break
            prev_pos = h.tell()
        h.seek(prev_pos)
        pos = h.tell()
        h.truncate(prev_pos + 1)
        h.seek(0, os.SEEK_END)
        assert h.tell() == pos

    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        simulation.open(".rxener", append=continue_sim)
        if continue_sim:
            self.__truncate(simulation)
        else:
            write_energies("", ".rxener", simulation, True)

    def __write_pos(self, title, sim: Simulation):
        if self.write_coords:
            sim.save_geometry(sim.request_path(f"_{title}{sim.current_step}{self.ext}"))

    def pre_modification(self, simulation: Simulation):
        # minimizations happen in "on_reaction", this is guaranteed to be before it
        write_energies("pre-reaction", ".rxener", simulation)
        self.__write_pos("premin", simulation)

    def post_reaction(self, simulation):
        write_energies("post-minimization", ".rxener", simulation)
        self.__write_pos("postmin", simulation)
