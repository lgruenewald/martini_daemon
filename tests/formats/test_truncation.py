from martini_daemon import Simulation, Reporter, FragCountReporter, VariablesReporter, ReactionReporter, \
    ReactionEnergyReporter, TrajectoryReporter
from typing import Type
import pytest
import os
import shutil

class FakeSimulation(Simulation):
    def __init__(self, base: str, reporters: list[Reporter], continue_sim: bool):
        # should keep current step inclusive
        self.current_step = 3000
        self.trajectory_frame = 4
        self.time_ps = 0.02 * 3000
        self.base = base
        self.reporters = reporters
        for r in reporters:
            r.on_simulation_start(self, continue_sim)

    def request_path(self, suffix: str, copy: bool = False) -> str:
        return self.base + suffix

    def finish(self) -> None:
        for r in self.reporters:
            r.on_simulation_finish(self)

def compare(a: str, b: str, mode: str) -> None:
    with open(a, mode) as f:
        a_content = f.read()
    with open(b, mode) as f:
        b_content = f.read()
    assert a_content.strip() == b_content.strip()

@pytest.fixture
def rootdir(request: pytest.FixtureRequest) -> str:
    """Fixture to get root directory."""
    return os.path.dirname(request.path)

formats = [
    (".frags", FragCountReporter, "r"),
    (".ener", VariablesReporter, "r"),
    (".reactions", ReactionReporter, "r"),
    (".rxener", ReactionEnergyReporter, "r"),
    (".xtc", TrajectoryReporter, "rb"),
]

@pytest.mark.parametrize("format_", formats)
def test_truncation(rootdir: str, format_: tuple[str, Type[Reporter]]) -> None:
    extension, reporter, mode = format_
    os.chdir(rootdir)
    os.chdir("to_truncate")
    shutil.copy(f"full{extension}", f"partial{extension}")

    sim = FakeSimulation("partial", [reporter()], continue_sim=True)
    sim.finish()
    compare(f"partial{extension}", f"reference{extension}", mode)
    os.remove(f"partial{extension}")
