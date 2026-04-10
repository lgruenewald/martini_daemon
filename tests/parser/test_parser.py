from typing import Any
import pytest
import os
from martini_daemon import Directive, Parser, TokenList


@pytest.fixture
def rootdir(request):
    return os.path.dirname(request.path)


class Root(Directive):
    def __init__(self):
        super().__init__(None, "", 0)
        self.root_lines = []
        self.a_s = []
        self.b_s = []
        self.finish_events = []

    def where(self) -> tuple[str, int]:
        raise NotImplementedError()

    def line(self, tokens: TokenList) -> None:
        self.root_lines.append([tokens.unwrap(i, "raw") for i in range(len(tokens))])

    def finish(self):
        pass

    @classmethod
    def is_mandatory(cls):
        raise NotImplementedError()

    @classmethod
    def is_unique(cls):
        raise NotImplementedError()

    @classmethod
    def is_valid_parent(cls, parent: Any) -> bool:
        raise NotImplementedError()

    @classmethod
    def get_name(cls) -> str:
        return "root"


class A(Directive):
    def __init__(self, parent: Any, path: str, line_num: int):
        super().__init__(parent, path, line_num)
        parent.a_s.append(self)
        self.lines = []

    def line(self, tokens: TokenList) -> None:
        self.lines.append([tokens.unwrap(i, "raw") for i in range(len(tokens))])

    def finish(self):
        self.parent.finish_events.append("a")

    @classmethod
    def is_mandatory(cls):
        return True

    @classmethod
    def is_unique(cls):
        return True

    @classmethod
    def is_valid_parent(cls, parent: Any) -> bool:
        return type(parent) is Root

    @classmethod
    def get_name(cls) -> str:
        return "a"


class B(Directive):

    def __init__(self, parent: Any, path: str, line_num: int):
        super().__init__(parent, path, line_num)
        parent.b_s.append(self)
        self.lines = []
        self.finish_events = parent.finish_events

    def line(self, tokens: TokenList) -> None:
        self.lines.append([tokens.unwrap(i, "raw") for i in range(len(tokens))])

    def finish(self):
        self.parent.finish_events.append("b")

    @classmethod
    def is_mandatory(cls):
        return False

    @classmethod
    def is_unique(cls):
        return False

    @classmethod
    def is_valid_parent(cls, parent: Any) -> bool:
        return type(parent) is Root

    @classmethod
    def get_name(cls) -> str:
        return "b"


class C(Directive):

    def __init__(self, parent: Any, path: str, line_num: int):
        super().__init__(parent, path, line_num)
        self.lines = parent.lines
        self.lines.append(["START OF C"])

    def line(self, tokens: TokenList) -> None:
        self.lines.append([tokens.unwrap(i, "float") for i in range(len(tokens))])

    def finish(self):
        self.parent.finish_events.append("c")

    @classmethod
    def is_mandatory(cls):
        return False

    @classmethod
    def is_unique(cls):
        return False

    @classmethod
    def is_valid_parent(cls, parent: Any) -> bool:
        return type(parent) is B

    @classmethod
    def get_name(cls) -> str:
        return "c"


def parse(path) -> tuple[bool, Root]:
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


def test_parser(rootdir: str):
    """
    Test Parser() with manually added directives.

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
def test_should_fail(x: str, rootdir):
    """
    Test the Directive API's guarantee requirements.
    """
    ok, root = parse(x)
    assert not ok
