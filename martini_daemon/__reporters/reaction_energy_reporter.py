from ..__reporter import Reporter
from ..__simulation import Simulation
from .variables_reporter import write_energies


class ReactionEnergyReporter(Reporter):
    def __init__(self, write_coords=False, ext=".gro"):
        """
        Reporter that will write energies before and after a reaction.

        :param write_coords: if set to True, it will print .gro files pre and post minimization
        :param ext: extension for writing the geometries. Set to ".xyz" if xyz files are desired.
        """
        self.write_coords = write_coords
        self.ext = ext

    def on_simulation_start(self, simulation: Simulation) -> None:
        simulation.open(".rxener")
        write_energies("", ".rxener", simulation, True)

    def __write_pos(self, title, sim: Simulation):
        if self.write_coords:
            sim.save_geometry(sim.request_path(f"_{title}{sim.current_step}{self.ext}"))

    def on_reaction(self, simulation: Simulation, _):
        write_energies("pre-reaction", ".rxener", simulation)
        self.__write_pos("premin", simulation)

    def post_reaction(self, simulation):
        write_energies("post-minimization", ".rxener", simulation)
        self.__write_pos("postmin", simulation)
