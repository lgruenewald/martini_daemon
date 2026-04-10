from ..__formats import TrajectoryWriter
from ..__reporter import Reporter
from ..__simulation import Simulation


class XTCReporter(Reporter):
    def __init__(self):
        self.path: str | None = None
        self.writer: TrajectoryWriter | None = None

    def on_simulation_start(self, simulation: Simulation) -> None:
        self.path = simulation.request_path(".xtc")
        self.writer = TrajectoryWriter(self.path)

    def on_trajectory_frame(self, simulation: Simulation) -> None:
        assert self.writer is not None
        assert simulation.__context is not None, (
            "Simulation context is None. Was it constructed with "
        )
        pos, box = simulation.__context.get_positions()
        self.writer.write_frame(simulation.current_step, simulation.time_ps, box, pos)

    def on_simulation_finish(self, simulation: Simulation) -> None:
        assert self.writer is not None
        self.writer.finish()
