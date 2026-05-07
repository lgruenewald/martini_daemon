import os
from typing import TextIO

from ..__simulation import Reporter, Simulation


def write_energies(
    title: str, handle: TextIO, sim: Simulation, first: bool = False
) -> None:
    if first:
        handle.write(
            "# Entry type,Simulation step,N,Kinetic energy (kJ/mol),Potential energy (kJ/mol),"
            "Total energy (kJ/mol),"
            "Temperature (Kelvin),"
            "Box X (nm),Box Y (nm),Box Z (nm),Volume (nm^3)\n",
        )
        return
    n = sim.system.num_atoms()
    ke, pe, te = sim.context.get_energies()
    degrees_of_freedom = sim.system.get_number_of_degrees_of_freedom()
    t = ke / degrees_of_freedom / 0.008314 * 2
    _, box = sim.context.get_positions()
    box_x = box.a[0]
    box_y = box.b[1]
    box_z = box.c[2]
    v = box_x * box_y * box_z
    handle.write(
        f"{title},{sim.current_step},{n},{ke},{pe},{te},{t},{box_x},{box_y},{box_z},{v}\n",
    )
    handle.flush()


def truncate_energies(handle: TextIO, sim: Simulation) -> None:
    handle.seek(0, os.SEEK_SET)
    prev_pos = 0
    while len(line := handle.readline()) > 0:
        if len(line) > 0 and line[0] != "#":
            frame = int(line.split(",")[1])
            if frame > sim.current_step:
                break
        prev_pos = handle.tell()
    handle.truncate(prev_pos)
    handle.seek(0, os.SEEK_END)


class VariablesReporter(Reporter):
    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        self.path = simulation.request_path(".ener", continue_sim=continue_sim)
        truncate = os.path.exists(self.path)
        self.handle = open(self.path, "r+" if truncate else "w")  # noqa: SIM115

        if truncate:
            truncate_energies(self.handle, simulation)
        else:
            self.handle.seek(0, os.SEEK_END)
            write_energies("", self.handle, simulation, True)

    def on_simulation_finish(self, simulation: Simulation) -> None:
        self.handle.close()

    def on_trajectory_frame(self, simulation: Simulation) -> None:
        write_energies("Trajectory frame", self.handle, simulation)
