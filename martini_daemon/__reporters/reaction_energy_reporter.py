import os

from ..__reporter import Reporter
from ..__simulation import Simulation
from .variables_reporter import truncate_energies, write_energies


class ReactionEnergyReporter(Reporter):
    def __init__(self, write_coords=False, ext=".gro"):
        """Reporter that will write energies before and after a reaction.

        :param write_coords: if set to True, it will print .gro files pre and post minimization
        :param ext: extension for writing the geometries. Set to ".xyz" if xyz files are desired.
        """
        self.write_coords = write_coords
        self.ext = ext

    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        self.path = simulation.request_path(".rxener", copy=continue_sim)
        truncate = continue_sim and os.path.exists(self.path)
        self.handle = open(self.path, "r+" if truncate else "w")  # noqa: SIM115

        if truncate:
            truncate_energies(self.handle, simulation)
        else:
            write_energies("", self.handle, simulation, True)

    def __write_pos(self, title, sim: Simulation):
        if self.write_coords:
            sim.save_geometry(sim.request_path(f"_{title}{sim.current_step}{self.ext}"))

    def pre_modification(self, simulation: Simulation):
        # minimizations happen in "on_reaction", this is guaranteed to be before it
        write_energies("pre-reaction", self.handle, simulation)
        self.__write_pos("premin", simulation)

    def post_reaction(self, simulation):
        write_energies("post-minimization", self.handle, simulation)
        self.__write_pos("postmin", simulation)

    def on_simulation_finish(self, simulation) -> None:
        self.handle.close()
