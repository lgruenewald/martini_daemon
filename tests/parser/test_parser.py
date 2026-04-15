"""Test the base Parsing of Martini Daemon, and the C preprocessor subset implemented."""

import os

import pytest

from martini_daemon import Directive, Parser, TokenList


@pytest.fixture
def rootdir(request: pytest.FixtureRequest) -> str:
    """Get the base directory where this test is located."""
    return os.path.dirname(request.path)


class Root(Directive):
    """Example root directive."""

    def __init__(self) -> None:
        """Create an example root directive."""
        super().__init__(None, "", 0)
        self.root_lines = []
        self.a_s = []
        self.b_s = []
        self.finish_events = []

    def where(self) -> tuple[str, int]:
        """Return where the directive started."""
        raise NotImplementedError()

    def line(self, tokens: TokenList) -> None:
        """Parse one line of input."""
        self.root_lines.append([tokens.unwrap(i, "raw") for i in range(len(tokens))])

    def finish(self) -> None:
        """Finish parsing this directive."""
        pass

    @classmethod
    def is_mandatory(cls) -> bool:
        """Get whether this directive is mandatory."""
        raise NotImplementedError()

    @classmethod
    def is_unique(cls) -> bool:
        """Get whether this directive is unique."""
        raise NotImplementedError()

    @classmethod
    def is_valid_parent(cls, parent: Directive) -> bool:
        """Get whether this object is a valid parent for this directive."""
        raise NotImplementedError()

    @classmethod
    def get_name(cls) -> str:
        """Get the name of this directive."""
        return "root"


class A(Directive):
    """Mandatory, unique example directive."""

    def __init__(self, parent: Directive, path: str, line_num: int) -> None:
        """Create an example child directive."""
        super().__init__(parent, path, line_num)
        assert isinstance(parent, Root)
        parent.a_s.append(self)
        self.lines = []

    def line(self, tokens: TokenList) -> None:
        """Parse one line of input."""
        self.lines.append([tokens.unwrap(i, "raw") for i in range(len(tokens))])

    def finish(self) -> None:
        """Finish parsing this directive."""
        self.parent.finish_events.append("a")

    @classmethod
    def is_mandatory(cls) -> bool:
        """Get whether this directive is mandatory."""
        return True

    @classmethod
    def is_unique(cls) -> bool:
        """Get whether this directive is unique."""
        return True

    @classmethod
    def is_valid_parent(cls, parent: Directive) -> bool:
        """Get whether this object is a valid parent for this directive."""
        return type(parent) is Root

    @classmethod
    def get_name(cls) -> str:
        """Get the name of this directive."""
        return "a"


class B(Directive):
    """Non-unique directive at Root."""

    def __init__(self, parent: Directive, path: str, line_num: int) -> None:
        """Create an example child directive."""
        super().__init__(parent, path, line_num)
        assert isinstance(parent, Root)
        parent.b_s.append(self)
        self.lines = []
        self.finish_events = parent.finish_events

    def line(self, tokens: TokenList) -> None:
        """Parse one line of input."""
        self.lines.append([tokens.unwrap(i, "raw") for i in range(len(tokens))])

    def finish(self) -> None:
        """Finish parsing this directive."""
        self.parent.finish_events.append("b")

    @classmethod
    def is_mandatory(cls) -> bool:
        """Get whether this directive is mandatory."""
        return False

    @classmethod
    def is_unique(cls) -> bool:
        """Get whether this directive is unique."""
        return False

    @classmethod
    def is_valid_parent(cls, parent: Directive) -> bool:
        """Get whether this object is a valid parent for this directive."""
        return type(parent) is Root

    @classmethod
    def get_name(cls) -> str:
        """Get the name of this directive."""
        return "b"


class C(Directive):
    """Child directive within B."""

    def __init__(self, parent: Directive, path: str, line_num: int) -> None:
        """Create a deeper example child directive of B."""
        super().__init__(parent, path, line_num)
        assert isinstance(parent, B)
        self.lines = parent.lines
        self.lines.append(["START OF C"])

    def line(self, tokens: TokenList) -> None:
        """Parse one line of input."""
        self.lines.append([tokens.unwrap(i, "float") for i in range(len(tokens))])

    def finish(self) -> None:
        """Finish parsing this directive."""
        self.parent.finish_events.append("c")

    @classmethod
    def is_mandatory(cls) -> bool:
        """Get whether this directive is mandatory."""
        return False

    @classmethod
    def is_unique(cls) -> bool:
        """Get whether this directive is unique."""
        return False

    @classmethod
    def is_valid_parent(cls, parent: Directive) -> bool:
        """Get whether this object is a valid parent for this directive."""
        return type(parent) is B

    @classmethod
    def get_name(cls) -> str:
        """Get the name of this directive."""
        return "c"


def parse(path: str) -> tuple[bool, Root]:
    """Parse an example root directive at path."""
    root = Root()
    parser = Parser(
        root=root,
        path=path,
        include_dirs=["includes"],
        defines={"TEST_DEFINE": "", "TEST_DEFINE2": "2"},
    )
    parser.add_directive(A)
    parser.add_directive(B)
    parser.add_directive(C)
    return parser.parse(), root


def test_parser(rootdir: str) -> None:
    """Test Parser() with manually added directives.

    Relevant file: test.ini.
    """
    os.chdir(rootdir)
    ok, root = parse("test.ini")
    assert ok
    root_lines = [
        ["root", "line1"],
        [
            "root",
            "line2",
            "token1",
            "token2",
            "token3",
            '"token4 with space"',
            "< token5 with spaces >",
        ],
    ]
    a_s = [
        [
            ["line1", "val1"],
            ["line2", "val2"],
        ]
    ]
    b_s = [
        [
            ["b1", "line"],
            ["START OF C"],
            [1.0, 2.0],
        ],
        [
            ["b2", "line"],
            ["b2", "line", "again"],
            ["START OF C"],
            [-1, -2],
            ["START OF C"],
            [3.5, 0.9],
        ],
    ]
    finish_events = ["c", "b", "c", "c", "b", "a"]
    for got, exp in zip(finish_events, root.finish_events):
        assert got == exp
    for got_list, exp_list in zip(root.root_lines, root_lines):
        for got, exp in zip(got_list, exp_list):
            assert got == exp
    for got_a, exp_list in zip(root.a_s, a_s):
        for got_line, exp_line in zip(got_a.lines, exp_list):
            for got, exp in zip(got_line, exp_line):
                assert got == exp
    for got_b, exp_list in zip(root.b_s, b_s):
        for got_line, exp_line in zip(got_b.lines, exp_list):
            for got, exp in zip(got_line, exp_line):
                assert got == exp


should_fail = ["fail_mandatory.ini", "fail_unique.ini", "fail_wrong_parent.ini"]


@pytest.mark.parametrize("x", should_fail)
def test_should_fail(x: str, rootdir: str) -> None:
    """Test the Directive API's guarantee requirements."""
    ok, root = parse(x)
    assert not ok
