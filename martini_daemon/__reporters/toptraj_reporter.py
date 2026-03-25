from ..__reporter import Reporter
from ..__simulation import Simulation
from ..__formats import TopTrajWriter
from ..__core import collect_bonds

class ToptrajReporter(Reporter):
    def __init__(self):
        self.writer: TopTrajWriter | None = None

    def on_simulation_start(self, simulation: Simulation):
        self.writer = TopTrajWriter(
            simulation.request_path(".toptraj"),
            simulation.system.additional_data.get("title"),
            simulation.system.initial_molecules
        )

    def on_trajectory_frame(self, simulation):
        self.writer.new_frame(
            simulation.trajectory_frame,
            simulation.current_step,
            simulation.time_ns,
            simulation.system.atom_count()
        )
        self.writer.register_frame_atoms(
            simulation.system.get_atom_names(),
            simulation.system.get_res_names(),
            simulation.system.get_res_ids(),
            simulation.system.get_types(),
            simulation.system.get_charges(),
            simulation.system.get_masses()
        )
        self.writer.register_frame_bonds(
            collect_bonds(simulation.system, ["vsite", "bond"]),
        )
        self.writer.write_frame()

    def on_simulation_finish(self, simulation):
        self.writer.finish()

