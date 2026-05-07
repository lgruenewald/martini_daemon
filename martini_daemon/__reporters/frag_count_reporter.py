import os

from ..__reporter import Reporter
from ..__simulation import Simulation


class FragCountReporter(Reporter):
    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        self.path = simulation.request_path(".frags", continue_sim=continue_sim)
        truncate = continue_sim and os.path.exists(self.path)
        self.handle = open(self.path, "r+" if truncate else "w")  # noqa: SIM115

        if truncate:
            truncate_to = simulation.current_step
            tell = self.handle.tell()
            while (line := self.handle.readline()) != "":
                c_step = int(line.split(",")[0].strip("Step:"))
                if c_step > truncate_to:
                    break
                tell = self.handle.tell()
            self.handle.truncate(tell)
            self.handle.seek(0, os.SEEK_END)

    def on_trajectory_frame(self, simulation: Simulation):
        self.handle.write(
            f"Step:{simulation.current_step},"
            + ",".join(
                [f"{k}:{v}" for k, v in simulation.top.frag_list.frag_counts.items()]
            )
            + "\n"
        )

    def on_simulation_finish(self, simulation) -> None:
        self.handle.close()

    def interactive_line(self, simulation) -> str:
        return f"fragments: {simulation.top.frag_list.num_fragments()}"

    def finish(self, simulation: Simulation) -> None:
        self.handle.close()
