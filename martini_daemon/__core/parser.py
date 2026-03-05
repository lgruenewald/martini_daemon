import re
from enum import Enum
import os
import traceback
from sys import stderr
from typing import Type
from .token import Token
from .token_list import TokenParseError, TokenList
from .directive import Directive


class ParseError(Exception):
    """
    Generic exception raised during parsing.
    """
    def __init__(self, message: str):
        self.message = message


class IfStackElem(Enum):
    DoBranch = 1
    SkipBranch = 2
    SkippedIf = 3
    Root = 4


class Parser:
    """
    Gromacs-style .ini format parser (used for .top/.itp files).

    .. note::
        Documentation below this line is kept as a reference for development.
        Most users do not need to interact with this API.

    A parser for gromacs topology file-like config files with limited C
    style preprocessing, that should be sufficient for most input files.
    You can add your own levels (Directive types). Your Directive types
    get instantiated when a level starts, with the parent directive instance
    as the constructor argument (often representing the System, or the currently
    edited Molecule). The line is method called on every line with TokenList as
    an argument (you can call TokenList.unwrap() to get individual tokens
    in your desired type). The finish() method is called when the directive ends.
    See :doc:`Directive</autoapi/martini_daemon/Directive>` class for details.

    All exceptions inside callbacks will be caught and re-raised,
    so that the specific line and file they occurred on can be printed to
    stdout, and the traceback suppressed to reduce clutter.
    """

    def __init__(self, root: Directive, path: str, include_dirs: list[str] | None = None, defines: dict[str, str] | None = None):
        """
        :param root: gets passed as the parent argument for directives with no parent directive. Can be useful to put the main root object being constructed here. Should be a directive instance.
        :param path: the path to the file to be parsed.
        :param include_dirs: list of directories to be searched if #include fails to find a file in the current dir.
        :param defines: dict[str, str] of keys and values for token replacements by the limited C preprocessor impl.
        """
        # stores types, not instances
        self.__directives: dict[str, Type[Directive]] = {}
        self.__line_num: int = 0
        self.__path: str = path
        self.__defines: dict[str, str] = defines or {}
        self.__include_dirs: list[str] = include_dirs or []
        self.__included: set[str] = set()
        self.__past_directives_at_root: set[str] = set()
        # stores instances
        # root is handled separately, as it's not guaranteed to be an instance of Directive
        self.__directive_stack: list[Directive] = [root]
        # state
        self.__done = False

    def add_directive(self, directive: Type[Directive]) -> None:
        """
        Add a new directive to this Parser.

        :param directive: the directive to add.
        """
        name = directive.get_name()
        assert name not in self.__directives, (
            f"Directive [{name}] was already defined."
        )
        self.__directives[name] = directive

    def parse(self) -> bool:
        """
        Performs the parsing and mutates the root argument passed during construction.

        :returns: Whether parsing succeeded. If not, the error message was already printed to stderr.

        Note: can only be called once per Parser instance.
        """
        if self.__done:
            print(
                "Parser.parse() must be called only once.",
                file=stderr
            )
            return False
        self.__done = True
        try:
            self.__parse(self.__path)
            for key, directive in self.__directives.items():
                if directive.is_mandatory() and key not in self.__past_directives_at_root:
                    raise ParseError(
                        f"Mandatory directive [{key}] not found."
                    )
            return True
        except TokenParseError as pe:
            self.__error_message_location(
                pe.message,
                pe.token.path, pe.token.line_num,
                pe.token.line,
                pe.token.start, pe.token.end
            )
            print()
        except ParseError as pe:
            self.__error_message_location(
                pe.message,
                self.__path, self.__line_num
            )
            print()
        except Exception:
            # yes it's broad, but we want to print where it happened
            traceback.print_exc()
            self.__error_message_location(
                "Exception occurred while parsing.",
                self.__path, self.__line_num
            )
        return False

    @staticmethod
    def __error_message_location(message: str, path: str, line_num: int, line=None, start=None, end=None) -> None:
        print("Parsing error:", file=stderr)
        print(message, file=stderr)
        print(f"In file {path} at line {line_num + 1}.", file=stderr)

        if line is not None and (start is None or end is None):
            print(line, file=stderr)
        else:
            print(
                f"{line[0:start]}"
                f"\033[1;33m{line[start:end]}\033[0m"
                f"{line[end:]}",
                file=stderr
            )

    # A token is either:
    # '['
    # ']'
    # "" enclosed any string
    # <> enclosed any string
    # contiguous array of characters that is not ASCII whitespace, [, ], ", < or >
    __token_pat = re.compile(r'\[|]|"[^"]*"|<[^>]*>|[^ \t\n\r\f\v\[\]"<>]+')

    def __tokenize(self, line: str) -> TokenList:
        """
        Splits a line-up into a list of tokens.

        * Assumes lines have been made whole already, and that comments were removed.
        * Separates based on whitespace.
        * Performs #define replacements.
        """

        tokens = [
            Token(
                content=line[match.start():match.end()],
                line=line,
                line_num=self.__line_num,
                path=self.__path,
                start=match.start(),
                end=match.end(),
            )
            for match in self.__token_pat.finditer(line)
        ]
        return TokenList(line, tokens, self.__defines)

    def __start_directive(self, dirname: str) -> None:
        directive_type = self.__directives.get(dirname)
        if directive_type is None:
            raise ParseError(
                f"Directive {dirname} not found.",
            )
        # == search for the right parent ==
        # should never pop root
        while len(self.__directive_stack) > 1:
            # find the outermost valid parent
            if directive_type.is_valid_parent(self.__directive_stack[-1]):
                break
            else:
                # directive ends previous directive -- call finishers
                # we need good error messages from finish hooks,
                # so temporarily edit path and line num to where the directive started
                old_ln = self.__line_num
                old_path = self.__path
                old_dir = self.__directive_stack.pop()
                self.__path, self.__line_num = old_dir.where()
                # finisher
                old_dir.finish()
                # since there's only one exception displayed at a time, if we get here
                # .finish() did not raise any exceptions
                self.__path, self.__line_num = old_path, old_ln
        parent = self.__directive_stack[-1]
        if not directive_type.is_valid_parent(parent):
            raise ParseError(
                f"Parent directive search for {dirname} failed. {parent.get_name()} is invalid parent."
            )
        # == check uniques (currently only for root) ==
        if len(self.__directive_stack) == 1:
            if directive_type.is_unique() and dirname in self.__past_directives_at_root:
                raise ParseError(
                    f"Unique directive {dirname} present more than once."
                )
            self.__past_directives_at_root.add(dirname)
        # == instantiate ==
        directive = directive_type(parent, self.__path, self.__line_num)
        self.__directive_stack.append(directive)

    def __parse(self, path: str) -> None:
        """
        Parses a single file at path.

        Recursively calls itself through #includes.
        """

        self.__included.add(path)
        old_path = self.__path
        self.__path = path

        # if stack is per file, but this is not a problem in most .top files
        if_stack = [IfStackElem.Root]

        with (open(path, "r") as f):
            # make lines whole through \
            cumulative = ""
            for i, line in enumerate(f):
                self.__line_num = i
                # ignore comments first to enable
                # blah blah \; comment
                # handle ignoring line endings
                if len(line) > 0 and line[-1] == "\\":
                    # must maintain whitespace, let's not strip
                    # note this will have implications for ; comment \
                    # but I think people shouldn't put \ after comment lines anyway
                    cumulative += line
                    continue
                if len(cumulative) > 0:
                    line = cumulative + line
                    cumulative = ""

                # remove comments and leading/trailing whitespace
                line = re.sub(r";.*", "", line).strip()
                # tokenize, doesn't resolve #define's yet
                token_list = self.__tokenize(line)

                # ignore empty lines
                # only after checking for \ -- means that empty lines
                # also need \ to keep continuing one long line
                if len(token_list) == 0:
                    continue

                # handle if/else logic before other things
                token0 = token_list[0]
                token1 = token_list[1] if len(token_list) > 1 else None
                if token0.content[0] == "#":
                    match token0.content:
                        case "#ifdef":
                            if token1 is None:
                                raise TokenParseError(
                                    token0, "#ifdef takes one argument."
                                )
                            if if_stack[-1] in {
                                IfStackElem.DoBranch,
                                IfStackElem.Root
                            }:
                                if token1.content in self.__defines:
                                    if_stack.append(IfStackElem.DoBranch)
                                else:
                                    if_stack.append(IfStackElem.SkipBranch)
                            else:
                                if_stack.append(IfStackElem.SkippedIf)
                            continue
                        case "#ifndef":
                            if token1 is None:
                                raise TokenParseError(
                                    token0, "#ifndef takes one argument."
                                )
                            if if_stack[-1] in {
                                IfStackElem.DoBranch,
                                IfStackElem.Root
                            }:
                                if token1.content in self.__defines:
                                    if_stack.append(IfStackElem.SkipBranch)
                                else:
                                    if_stack.append(IfStackElem.DoBranch)
                            else:
                                if_stack.append(IfStackElem.SkippedIf)
                            continue
                        case "#else":
                            if token1 is not None:
                                raise TokenParseError(
                                    token0, "#else takes no argument."
                                )
                            match if_stack[-1]:
                                case IfStackElem.DoBranch:
                                    if_stack[-1] = IfStackElem.SkipBranch
                                case IfStackElem.SkipBranch:
                                    if_stack[-1] = IfStackElem.DoBranch
                                case IfStackElem.SkippedIf:
                                    pass
                                case IfStackElem.Root:
                                    raise TokenParseError(
                                        token0, "#else unmatched."
                                    )
                            continue
                        case "#endif":
                            if token1 is not None:
                                raise TokenParseError(
                                    token0, "#endif takes no argument."
                                )
                            if if_stack[-1] == IfStackElem.Root:
                                raise TokenParseError(
                                    token0, "#endif unmatched."
                                )
                            if_stack.pop()
                            continue
                        case "#end":
                            # possibly common mistake?
                            raise TokenParseError(
                                token0, "Please use #endif."
                            )
                        case "#if" | "#elif":
                            raise TokenParseError(
                                token0, "Only #ifdef is implemented. #if, #elif is not."
                            )
                        case "#include" | "#define" | "#undef":
                            # preprocessor directives that are conditional
                            pass
                        case _:
                            # error at unknown macros
                            raise TokenParseError(
                                token0, "Unknown preprocessor directive."
                            )

                # everything below this only happens if the #ifdef/#else
                # says it should happen
                if if_stack[-1] in {
                    IfStackElem.SkipBranch,
                    IfStackElem.SkippedIf
                }:
                    continue

                if token0.content == "[":
                    # directives
                    if token_list[-1].content != "]":
                        raise TokenParseError(
                            token_list[-1], "Invalid directive, not closed by ']'."
                        )
                    if len(token_list) != 3:
                        raise TokenParseError(
                            token_list[-1],
                            "Invalid directive, three tokens expected: '[', directive name, ']'."
                        )

                    self.__start_directive(token_list.unwrap(1, "word"))
                elif token0.content[0] == "#":
                    match token0.content:
                        case "#include":
                            name = token_list.unwrap(
                                1,
                                "string",
                                error_msg="#include argument should be inside quotation marks or <>."
                            )
                            search_dirs = [
                                os.path.dirname(path),
                            ] + self.__include_dirs

                            found = False
                            for cdir in search_dirs:
                                new_path = os.path.join(cdir, name)
                                if os.path.isfile(new_path):
                                    if new_path in self.__included:
                                        raise TokenParseError(
                                            token1,
                                            f"Double inclusion of {new_path}."
                                        )
                                    self.__parse(new_path)
                                    found = True
                                    break
                            if not found:
                                raise TokenParseError(
                                    token1,
                                    f"File not found: {name}."
                                )
                        case "#define":
                            if len(token_list) not in {2, 3}:
                                raise TokenParseError(
                                    token0,
                                    "#define takes one or two tokens arguments. "
                                    "Note: Only single token -> single token mappings are supported. "
                                    "Preprocessor macros are not supported."
                                )
                            # let's be sane and not allow #define "blah.itp" "blah2.itp" and similar
                            key = token_list.unwrap(1, "word")
                            # now we can actually use an empty line, as long as it is in the dict #ifdef supports it
                            value = token_list.unwrap(2, "raw", default="")
                            self.__defines[key] = value
                        case "#undef":
                            if len(token_list) != 2:
                                raise TokenParseError(
                                    token0,
                                    "#undef takes one argument."
                                )
                            key = token_list.unwrap(1, "word")
                            if self.__defines.get(key):
                                self.__defines.pop(key)
                        case _:
                            # unreachable
                            assert False
                else:
                    # data lines
                    self.__directive_stack[-1].line(token_list)

        if len(if_stack) > 1:
            raise ParseError(
                "Unmatched #ifdef or #ifndef. #ifdef/#ifndef crossing file boundaries are not supported."
            )
        self.__path = old_path
