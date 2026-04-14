from ..__formats import TrajectoryWriter
from ..__reporter import Reporter
from ..__simulation import Simulation


class TrajectoryReporter(Reporter):
    def __init__(self, format_: str = ".xtc", backend: str | None = None) -> None:
        self.__format = format_
        self.__backend = backend
        self.path: str | None = None
        self.writer: TrajectoryWriter | None = None

    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        # TODO truncate
        if continue_sim:
            raise NotImplementedError()

        self.path = simulation.request_path(self.__format)
        self.writer = TrajectoryWriter(self.path, backend=self.__backend)

    def on_trajectory_frame(self, simulation: Simulation) -> None:
        assert self.writer is not None
        assert simulation.context is not None, (
            "Simulation context is None. Was it constructed with "
        )
        pos, box = simulation.context.get_positions()
        self.writer.write_frame(simulation.current_step, simulation.time_ps, box, pos)

    def on_simulation_finish(self, simulation: Simulation) -> None:
        assert self.writer is not None
        self.writer.finish()
