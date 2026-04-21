from ..__formats import TopTrajWriter
from ..__reporter import Reporter
from ..__simulation import Simulation


class TopTrajReporter(Reporter):
    def __init__(self):
        self.writer: TopTrajWriter | None = None

    def on_simulation_start(self, simulation: Simulation, continue_sim: bool):
        if continue_sim:
            raise NotImplementedError()  # TODO
        title = simulation.system.additional_data.get("title")
        assert type(title) is str
        self.writer = TopTrajWriter(
            simulation.request_path(".toptraj"),
            title,
            simulation.system.initial_molecules,
        )

    def on_trajectory_frame(self, simulation):
        assert self.writer is not None
        self.writer.new_frame(
            simulation.trajectory_frame,
            simulation.current_step,
            simulation.time_ps,
            simulation.system.num_atoms(),
        )
        self.writer.write_frame_atoms(
            simulation.system.get_atom_names(),
            simulation.system.get_res_names(),
            simulation.system.get_res_ids(),
            simulation.system.get_types(),
            simulation.system.get_charges(),
            simulation.system.get_masses(),
        )
        self.writer.write_frame_bonds(
            simulation.system.collect_bonds(["vsite", "bond"]),
        )
        self.writer.write_frame()

    def on_simulation_finish(self, simulation):
        assert self.writer is not None
        self.writer.finish()
