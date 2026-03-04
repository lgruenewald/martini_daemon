
"""
there needs to be some global state --> GromacsTopFormat should collect global state

this is done via decorators:

@directive
@force

these decorators should be added to classes
these classes should have predefined members based on what it is that gets marked
the user then does not ever instantiate these themselves, the simulation object instantiates them
the simulation object gets passed in __init__, exposing the same familiar API
@directive -- designed in a way where each instance of the directive in source files instantiates the class once
__init__ is called when the directive gets opened, line() gets called on a line of source, __del__ gets called
when the directive gets closed
@bond, @angle, @dihedral, @exclusion... - predefined directives that can be extended with custom forces,
since these are in __core, they can find the instance of the force in system to write to. the parsing logic for a single
force type can be implemented inside them
"""

from .__directive import Directive

class GromacsTopFormat:
    """
    Metadata about the .top format, allowing for the construction of parsers and generators of said format.
    """

    def __init__(self):
        """Create a new TopFormat parser."""
        self.__directives: list[Directive] = []

    def add_directive(self, directive: Directive) -> None:
        """Add a directive to TopFormat."""
        self.__directives.append(directive)
