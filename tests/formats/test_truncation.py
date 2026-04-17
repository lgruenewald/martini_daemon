from martini_daemon import Simulation, Reporter, FragCountReporter, VariablesReporter, ReactionReporter, ReactionEnergyReporter
from typing import Type
import pytest
import os
import shutil

class FakeSimulation(Simulation):
    def __init__(self, base: str, reporters: list[Reporter], continue_sim: bool):
        # should keep current step inclusive
        self.current_step = 3000
        self.base = base
        self.reporters = reporters
        for r in reporters:
            r.on_simulation_start(self, continue_sim)

    def request_path(self, suffix: str, copy: bool = False) -> str:
        return self.base + suffix

def compare(a: str, b: str) -> None:
    with open(a, "r") as f:
        a_content = f.read()
    with open(b, "r") as f:
        b_content = f.read()
    assert a_content.strip() == b_content.strip()

@pytest.fixture
def rootdir(request: pytest.FixtureRequest) -> str:
    """Fixture to get root directory."""
    return os.path.dirname(request.path)

formats = [
    (".frags", FragCountReporter),
    (".ener", VariablesReporter),
    (".reactions", ReactionReporter),
    (".rxener", ReactionEnergyReporter)
]

@pytest.mark.parametrize("format_", formats)
def test_truncation(rootdir: str, format_: tuple[str, Type[Reporter]]) -> None:
    extension, reporter = format_
    os.chdir(rootdir)
    os.chdir("to_truncate")
    shutil.copy(f"full{extension}", f"partial{extension}")

    FakeSimulation("partial", [reporter()], continue_sim=True)
    compare(f"partial{extension}", f"reference{extension}")
    os.remove(f"partial{extension}")
