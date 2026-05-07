"""Test the truncation of output files by comparing reference with truncated output."""

import os
import shutil

import pytest

from martini_daemon import (
    FragCountReporter,
    ReactionEnergyReporter,
    ReactionReporter,
    Reporter,
    Simulation,
    System,
    TopTrajReporter,
    TrajectoryReporter,
    VariablesReporter,
)


class FakeSystem(System):
    """Fake system class for testing."""

    def __init__(self) -> None:
        """Fake system object."""
        super().__init__()
        self.additional_data["title"] = "placeholder"

class FakeSimulation(Simulation):
    """Fake simulation class for testing."""

    def __init__(self, base: str, reporters: list[Reporter], continue_sim: bool) -> None:
        """Create fake simulation."""
        # should keep current step inclusive
        self.current_step = 3000
        self.trajectory_frame = 4
        self.time_ps = 0.02 * 3000
        self.base = base
        self.reporters = reporters
        self.system = FakeSystem()
        for r in reporters:
            r.on_simulation_start(self, continue_sim)

    def request_path(self, suffix: str, continue_sim: bool = False) -> str:
        """Request path to output file with suffix."""
        return self.base + suffix

    def finish(self) -> None:
        """Finish the fake simulation."""
        for r in self.reporters:
            r.on_simulation_finish(self)


def compare(a: str, b: str, mode: str) -> None:
    """Compare two files, opened with mode <mode>."""
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
    (".toptraj", TopTrajReporter, "rb"),
]


@pytest.mark.parametrize("format_", formats)
def test_truncation(rootdir: str, format_: tuple[str, type[Reporter], str]) -> None:
    """Truncate a source file, compare with reference."""
    extension, reporter, mode = format_
    os.chdir(rootdir)
    os.chdir("to_truncate")
    shutil.copy(f"full{extension}", f"partial{extension}")

    sim = FakeSimulation("partial", [reporter()], continue_sim=True)
    sim.finish()
    compare(f"partial{extension}", f"reference{extension}", mode)
    os.remove(f"partial{extension}")
