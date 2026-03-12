import pytest
import os
from typing import Any
from martini_daemon import GromacsTopFile, Directive, register_directive, TokenList

@pytest.fixture
def rootdir(request):
    return os.path.dirname(request.path)

@register_directive
class Custom(Directive):
    def __init__(self, parent: GromacsTopFile, path, num):
        super().__init__(parent, path, num)
        self.lines = []
        # not the cleanest code, since this test was written before GromacsTopFile's API was finished
        parent.injected = self

    def line(self, tokens: TokenList) -> None:
        self.lines.append(" ".join([
            tokens.unwrap(i, "raw") for i in range(len(tokens))
        ]))

    def finish(self):
        self.lines.append("Finished")

    @staticmethod
    def is_mandatory():
        return True

    @staticmethod
    def is_unique():
        return True

    @staticmethod
    def is_valid_parent(parent: Any) -> bool:
        return type(parent) is GromacsTopFile

    @staticmethod
    def get_name() -> str:
        return "custom"


def test_custom_directives(rootdir):
    """
    Test GromacsTopParser's @directive registration with custom directives.

    Relevant file: test_custom_directive.ini.
    """
    os.chdir(rootdir)
    top = GromacsTopFile("test_custom_directive.ini")
    expected = [
        "Will be parsed by custom parser", "Finished"
    ]
    for exp, got in zip(expected, top.injected.lines):
        assert exp == got