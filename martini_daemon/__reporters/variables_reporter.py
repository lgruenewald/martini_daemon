from ..__reporter import Reporter
from ..__simulation import Simulation


def write_energies(title, suffix, sim: Simulation, first=False):
    if first:
        sim.print(
            suffix,
            "# Entry type,Simulation step,N,Kinetic energy (kJ/mol),Potential energy (kJ/mol),"
            "Total energy (kJ/mol),"
            "Temperature (Kelvin),"
            "Box X (nm),Box Y (nm),Box Z (nm),Volume (nm^3)",
        )
        return
    n = sim.system.atom_count()
    ke, pe, te = sim.get_context().get_energies()
    degrees_of_freedom = sim.system.get_number_of_degrees_of_freedom()
    t = ke / degrees_of_freedom / 0.008314 * 2
    _, box = sim.get_context().get_positions()
    box_x = box.a[0]
    box_y = box.b[1]
    box_z = box.c[2]
    v = box_x * box_y * box_z
    sim.print(
        suffix,
        f"{title},{sim.current_step},{n},{ke},{pe},{te},{t},{box_x},{box_y},{box_z},{v}",
    )


class VariablesReporter(Reporter):
    def on_simulation_start(self, simulation: Simulation):
        simulation.open(".ener")
        write_energies("", ".ener", simulation, True)

    def on_trajectory_frame(self, simulation: Simulation):
        write_energies("Trajectory frame", ".ener", simulation)
