from ..__reporter import Reporter
from ..__simulation import Simulation
from ..__formats import TrajectoryWriter

class XTCReporter(Reporter):
    def __init__(self):
        self.path: str | None  = None
        self.writer: TrajectoryWriter | None = None

    def on_simulation_start(self, simulation):
        self.path = simulation.request_path(".xtc")
        self.writer = TrajectoryWriter(self.path)

    def on_trajectory_frame(self, simulation):
        pos, box = simulation.context.get_positions()
        self.writer.write_frame(simulation.current_step, simulation.time_ns, box, pos)

    def on_simulation_finish(self, simulation):
        self.writer.finish()