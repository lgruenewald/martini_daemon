from ..__formats import TopTrajWriter
from ..__simulation import Reporter, Simulation


class TopTrajReporter(Reporter):
    def __init__(self) -> None:
        """
        Create a TopTrajReporter.

        Will write a "topology trajectory" file to a file ending with `.toptraj`, written at each trajectory
        frame of the simulation.

        See :doc:`/toptraj`, :doc:`/autoapi/martini_daemon/TopTrajWriter`
        and :doc:`/autoapi/martini_daemon/TopTrajReader` for more information about the format.
        """
        self.writer: TopTrajWriter | None = None

    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        title = simulation.system.additional_data.get("title")
        assert type(title) is str
        init_mols = []
        for name, n in simulation.system.initial_molecules:
            n_atoms_per = len(simulation.system.molecule_types[name].atoms)
            init_mols.append((name, n, n_atoms_per))
        self.writer = TopTrajWriter(
            simulation.request_path(".toptraj", continue_sim=continue_sim),
            title,
            init_mols,
            simulation.system.get_res_names(),
            simulation.system.get_res_ids(),
            append=continue_sim,
            truncate=simulation.current_step,
        )

    def on_trajectory_frame(self, simulation: Simulation) -> None:
        assert self.writer is not None
        self.writer.new_frame(
            simulation.trajectory_frame,
            simulation.current_step,
            simulation.time_ps,
            simulation.system.num_atoms(),
        )
        self.writer.write_frame_atoms(
            simulation.system.get_atom_names(),
            simulation.system.get_types(),
            simulation.system.get_charges(),
            simulation.system.get_masses(),
        )
        self.writer.write_frame_bonds(
            simulation.system.collect_bonds(["vsite", "bond"]),
        )
        self.writer.write_frame()

    def on_simulation_finish(self, simulation: Simulation) -> None:
        assert self.writer is not None
        self.writer.finish()
