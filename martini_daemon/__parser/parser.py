import os
import re
import traceback
from enum import Enum
from sys import stderr

from ..__rust import Token, TokenList, Lexer, TokenParseException
from .directive import Directive


class ParseException(Exception):
    def __init__(self, message: str) -> None:
        """Create a generic exception raised during parsing."""
        self.message = message


class DirectiveException(Exception):
    def __init__(
        self, message: str, path: str, start_line: int, end_line: int | None
    ) -> None:
        """
        Create a Directive Exception.

        This exception is raised during the finish() method of directives.
        ParseExceptions get re-raised as this type automatically.
        """
        self.path = path
        self.end_line = end_line
        self.start_line = start_line
        self.message = message


class IfStackElem(Enum):
    DoBranch = 1
    SkipBranch = 2
    SkippedIf = 3
    Root = 4


class Parser:
    def __init__(
        self,
        root: Directive,
        path: str,
        include_dirs: list[str] | None = None,
        defines: dict[str, str] | None = None,
    ) -> None:
        """
        Create a generic Gromacs-style .ini format parser.

        See :doc:`/autoapi/martini_daemon/GromacsTopFile` for the parser specifically for GROMACS ``.top`` files,
        which inherits this class. For extending Martini Daemon with custom directives, also see :doc:`/extending`.

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

        :param root: gets passed as the parent argument for directives with no parent directive.
            Can be useful to put the main root object being constructed here. Should be a directive instance.
        :param path: the path to the file to be parsed.
        :param include_dirs: list of directories to be searched if #include fails to find a file in the current dir.
        :param defines: dict[str, str] of keys and values for token replacements by the limited C preprocessor impl.
        """
        # stores types, not instances
        self.__directives: dict[str, type[Directive]] = {}
        self.__path: str = path
        self.__lexer: Lexer = Lexer()
        self.__line_num: int = 0
        if defines is not None:
            for k, v in defines.items():
                self.__lexer.add_define(k, v)
        self.__include_dirs: list[str] = include_dirs or []
        self.__included: set[str] = set()
        self.__past_directives_at_root: set[str] = set()
        # stores instances
        # root is handled separately, as it's not guaranteed to be an instance of Directive
        self.__directive_stack: list[Directive] = [root]
        # state
        self.__done = False

    def add_directive(self, directive: type[Directive]) -> None:
        """
        Add a new directive to this Parser.

        :param directive: the directive to add.
        """
        name = directive.get_name()
        assert name not in self.__directives, f"Directive [{name}] was already defined."
        self.__directives[name] = directive
        for alias in directive.aliases:
            assert alias not in self.__directives, (
                f"Directive [{alias}] was already defined."
            )
            self.__directives[alias] = directive

    def parse(self) -> bool:
        """
        Perform the parsing and mutates the root argument passed during construction.

        :returns: Whether parsing succeeded. If not, the error message was already printed to stderr.

        Note: can only be called once per Parser instance.
        """
        if self.__done:
            print("Parser.parse() must be called only once.", file=stderr)
            return False
        self.__done = True
        try:
            self.__parse(self.__path)

            while len(self.__directive_stack) > 0:
                self.__directive_stack.pop().finish()

            for key, directive in self.__directives.items():
                if (
                    directive.is_mandatory()
                    and key not in self.__past_directives_at_root
                ):
                    raise ParseException(f"Mandatory directive [{key}] not found.")
            return True
        except TokenParseException as pe:
            self.__error_message_location(
                pe.message,
                pe.token.path,
                pe.token.line_num,
                pe.token.start,
                pe.token.end,
            )
        except ParseException as pe:
            self.__error_message_location(
                pe.message, self.__path, self.__line_num, None
            )
        except DirectiveException as de:
            self.__error_directive(
                de.message,
                de.path,
                de.start_line,
                de.end_line,
            )
        except Exception:
            # yes it's broad, but we want to print where it happened
            traceback.print_exc()
            self.__error_message_location(
                "Exception occurred while parsing.", self.__path, self.__line_num, None
            )
        return False

    @classmethod
    def __error_message_location(
        cls,
        message: str,
        path: str,
        line_num: int,
        start: int | None = None,
        end: int | None = None,
    ) -> None:
        """
        Print an error message.

        Three modes available:

        - ``line``, ``start``, ``end`` are all None - will print the line based on the file on disk
        - only ``line`` is not None - will print ``line`` as line content
        - ``line``, ``start``, ``end`` are all not None - will print ``line``, with the range ``start``:``end`` yellow
        """
        print(f"\033[1;33m{message}\033[0m", file=stderr)
        print(f"In file {path} at line {line_num + 1}.", file=stderr)

        with open(path) as file:
            lines = file.read().splitlines()
            if len(lines) <= line_num:
                return
            line = lines[line_num]

        if start is None or end is None:
            print(f"{line_num + 1}: \033[1;33m{line}\033[0m", file=stderr)
        else:
            assert line is not None
            print(
                f"{line_num + 1}: {line[0:start]}"
                f"\033[1;33m{line[start:end]}\033[0m"
                f"{line[end:]}",
                file=stderr,
            )

    @classmethod
    def __error_directive(
        cls, message: str, path: str, start: int, end: int | None
    ) -> None:
        """
        Print an error message.

        If end is not None, it will read the file and highlight the whole directive's
        text in yellow, with line numbers.
        """
        assert end is None or end > start

        print(f"\033[1;33m{message}\033[0m", file=stderr)
        print(f"In file {path} at line {start + 1}.", file=stderr)

        if end is not None:
            with open(path) as file:
                lines = file.read().splitlines()
                if len(lines) <= end:
                    return
                print(
                    "\n".join(
                        f"{i + 1}: \033[1;33m{line}\033[0m"
                        for i, line in zip(range(start, end), lines[start:end])
                    ),
                    file=stderr,
                )

    # A token is either:
    # '['
    # ']'
    # "" enclosed any string
    # <> enclosed any string
    # contiguous array of characters that is not ASCII whitespace, [, ], ", < or >
    __token_pat = re.compile(r'\[|]|"[^"]*"|<[^>]*>|[^ \t\n\r\f\v\[\]"<>]+')

    def __start_directive(self, dirname: str) -> None:
        directive_type = self.__directives.get(dirname)
        if directive_type is None:
            raise ParseException(
                f"Directive {dirname} not found.",
            )
        # == search for the right parent ==
        # should never pop root
        while len(self.__directive_stack) > 1:
            # find the outermost valid parent
            if directive_type.is_valid_parent(self.__directive_stack[-1]):
                break
            # directive ends previous directive -- call finishers
            old_dir = self.__directive_stack.pop()

            # info for potential exceptions
            directive_end = self.__line_num
            current_path = self.__path
            directive_path, directive_start = old_dir.where()
            # finisher
            try:
                old_dir.finish()
            except ParseException as pe:
                # we need good error messages from finish hooks, so let's re-raise ParseExceptions
                raise DirectiveException(
                    pe.message,
                    directive_path,
                    directive_start,
                    directive_end if current_path == directive_path else None,
                )
        parent = self.__directive_stack[-1]
        if not directive_type.is_valid_parent(parent):
            raise ParseException(
                f"Parent directive search for {dirname} failed. {parent.get_name()} is invalid parent."
            )
        # == check uniques (currently only for root) ==
        if len(self.__directive_stack) == 1:
            if directive_type.is_unique() and dirname in self.__past_directives_at_root:
                raise ParseException(
                    f"Unique directive {dirname} present more than once."
                )
            self.__past_directives_at_root.add(dirname)
        # == instantiate ==
        directive = directive_type(parent, self.__path, self.__line_num)
        self.__directive_stack.append(directive)

    def __parse(self, path: str) -> None:
        """
        Parse a single file at path.

        Recursively calls itself through #includes.
        """
        self.__included.add(path)
        old_path = self.__path
        self.__path = path

        # if stack is per file, but this is not a problem in most .top files
        if_stack = [IfStackElem.Root]

        lines: list[TokenList] = self.__lexer.tokenize_file(
            path
        )

        for token_list in lines:
            self.__line_num = token_list.line_num
            if len(token_list) == 0:
                continue
            # handle if/else logic before other things
            token0 = token_list[0]
            token1 = token_list[1] if len(token_list) > 1 else None
            if token0.content[0] == "#":
                match token0.content:
                    case "#ifdef":
                        if token1 is None or len(token_list) > 2:
                            raise TokenParseException(
                                token0, "#ifdef takes one argument."
                            )
                        if if_stack[-1] in {IfStackElem.DoBranch, IfStackElem.Root}:
                            if self.__lexer.get_define(token1.content) is not None:
                                if_stack.append(IfStackElem.DoBranch)
                            else:
                                if_stack.append(IfStackElem.SkipBranch)
                        else:
                            if_stack.append(IfStackElem.SkippedIf)
                        continue
                    case "#ifndef":
                        if token1 is None or len(token_list) > 2:
                            raise TokenParseException(
                                token0, "#ifndef takes one argument."
                            )
                        if if_stack[-1] in {IfStackElem.DoBranch, IfStackElem.Root}:
                            if self.__lexer.get_define(token1.content) is not None:
                                if_stack.append(IfStackElem.SkipBranch)
                            else:
                                if_stack.append(IfStackElem.DoBranch)
                        else:
                            if_stack.append(IfStackElem.SkippedIf)
                        continue
                    case "#else":
                        if token1 is not None:
                            raise TokenParseException(
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
                                raise TokenParseException(
                                    token0, "#else unmatched."
                                )
                        continue
                    case "#endif":
                        if token1 is not None:
                            raise TokenParseException(
                                token0, "#endif takes no argument."
                            )
                        if if_stack[-1] == IfStackElem.Root:
                            raise TokenParseException(token0, "#endif unmatched.")
                        if_stack.pop()
                        continue
                    case "#end":
                        # possibly common mistake?
                        raise TokenParseException(token0, "Please use #endif.")
                    case "#if" | "#elif":
                        raise TokenParseException(
                            token0, "Only #ifdef is implemented. #if, #elif is not."
                        )
                    case "#include" | "#define" | "#undef":
                        # preprocessor directives that are conditional
                        pass
                    case _:
                        # error at unknown macros
                        raise TokenParseException(
                            token0, "Unknown preprocessor directive."
                        )

            # everything below this only happens if the #ifdef/#else
            # says it should happen
            if if_stack[-1] in {IfStackElem.SkipBranch, IfStackElem.SkippedIf}:
                continue

            if token0.content == "[":
                # directives
                if token_list[-1].content != "]":
                    raise TokenParseException(
                        token_list[-1], "Invalid directive, not closed by ']'."
                    )
                if len(token_list) != 3:
                    raise TokenParseException(
                        token_list[-1],
                        "Invalid directive, three tokens expected: '[', directive name, ']'.",
                    )

                self.__start_directive(token_list.unwrap(1, "word"))
            elif token0.content[0] == "#":
                match token0.content:
                    case "#include":
                        # #include actually supports preprocessor replacements, unlike #ifdef or #define
                        self.__lexer.do_replacements(token_list)
                        if len(token_list) != 2:
                            raise TokenParseException(
                                token0, "#include takes one argument."
                            )
                        name = token_list.unwrap(
                            1,
                            "string",
                        )
                        search_dirs = [
                            os.path.dirname(path),
                        ] + self.__include_dirs

                        found = False
                        for cdir in search_dirs:
                            new_path = os.path.join(cdir, name)
                            if os.path.isfile(new_path):
                                if new_path in self.__included:
                                    raise TokenParseException(
                                        token1, f"Double inclusion of {new_path}."
                                    )
                                self.__parse(new_path)
                                found = True
                                break
                        if not found:
                            raise TokenParseException(
                                token1, f"File not found: {name}."
                            )
                    case "#define":
                        if token1 is None or len(token_list) < 2:
                            raise TokenParseException(
                                token0,
                                "#define takes an argument."
                            )
                        key = token1.content
                        value = " ".join([token_list[i].content for i in range(2, len(token_list))])
                        self.__lexer.add_define(key, value)
                    case "#undef":
                        if token1 is None or len(token_list) != 2:
                            raise TokenParseException(
                                token0, "#undef takes one argument."
                            )
                        key = token1.content
                        self.__lexer.rm_define(key)
                    case _:
                        # unreachable
                        assert False
            else:
                # data lines
                self.__lexer.do_replacements(token_list)
                # #defines could have made this line empty
                if len(token_list) > 0:
                    self.__directive_stack[-1].line(token_list)

        if len(if_stack) > 1:
            raise ParseException(
                "Unmatched #ifdef or #ifndef. #ifdef/#ifndef crossing file boundaries are not supported."
            )
        self.__path = old_path
