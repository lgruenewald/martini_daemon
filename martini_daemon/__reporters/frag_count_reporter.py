from ..__reporter import Reporter
from ..__simulation import Simulation


class FragCountReporter(Reporter):
    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        simulation.open(".frags", append=continue_sim)
        # TODO truncate

    def on_trajectory_frame(self, simulation: Simulation):
        simulation.print(
            ".frags",
            f"Step:{simulation.current_step},"
            + ",".join(
                [f"{k}:{v}" for k, v in simulation.top.frag_list.frag_counts.items()]
            ),
        )

    def interactive_line(self, simulation) -> str:
        return f"fragments: {simulation.top.frag_list.num_fragments()}"
