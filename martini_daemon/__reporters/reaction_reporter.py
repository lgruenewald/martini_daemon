import os
import re

from ..__reporter import Reporter
from ..__rust import Fragment
from ..__simulation import Simulation


class ReactionReporter(Reporter):
    """A reporter that reports all reactions to <name>.reactions"""

    def __init__(self):
        """A reporter that reports all reactions to <name>.reactions.

        The created file has a text format, where every line is a reaction.
        First, the frame number and reaction name are separated by a comma,
        then, each reactant is separated by a semicolon. Each reactant
        will have its frag name, internal frag id, and atom indices
        (-1 for missing optional or forbidden atoms) printed.

        The residue IDs involved are also included as a comma separated list, within parentheses, separately
        for each reactant.

        Example:
        frame,reaction_name;reactant1_name,reactant1_id(res:resid1,...residn),atoms...;...reactantn_name(res:resid1,...residn),reactantn_id,atoms...

        """
        self.reactions = 0

    def __truncate(self, sim: Simulation):
        h = sim.get_handle(".reactions")
        h.seek(0, os.SEEK_SET)
        prev_pos = 0
        while (line := h.readline()) != b"":
            frame = int(line.decode("utf-8").split(",")[0])
            if frame > sim.current_step:
                break
            prev_pos = h.tell()
        h.truncate(prev_pos)
        h.seek(0, os.SEEK_END)

    def on_simulation_start(self, simulation, continue_sim: bool):
        self.handle = open(
            simulation.request_path(".reactions", copy=continue_sim), "r+"
        )
        if continue_sim:
            self.__truncate(simulation)
        else:
            self.handle.seek(0, os.SEEK_END)
            simulation.print(
                ".reactions",
                "# sim step,reaction_name;"
                + "reactant1_name,reactant1_id(res:resid1,...residn),atoms...;..."
                + "reactantn_name,reactantn_id(res:resid1,...residn),atoms...;",
            )

    def on_simulation_finish(self, simulation) -> None:
        self.handle.close()

    @staticmethod
    def __get_resids(sim: Simulation, atoms):
        res = set()
        for atom in atoms:
            if atom != -1:
                res.add(f"{sim.system.get_res_ids()[atom]}")
        return "(res:" + ",".join(res) + ")"

    # need to be post, as the modification algorithm can reject some reactions
    def on_reaction(
        self, simulation: Simulation, reactions: list[tuple[str, list[Fragment]]]
    ):
        for rx, frags in reactions:
            simulation.print(
                ".reactions",
                f"{simulation.current_step},{rx};"
                + ";".join(
                    [
                        f"{frag.name},{frag.frag_id}"
                        f"{self.__get_resids(simulation, frag.atoms)},"
                        + ",".join([f"{atom}" for atom in frag.atoms])
                        for frag in frags
                    ]
                ),
            )
        self.reactions += len(reactions)

    def interactive_line(self, simulation) -> str:
        return f"reactions: {self.reactions}"

    @staticmethod
    def read_reactions(path) -> list[tuple[int, str, list[list[int]]]]:
        """.reactions format reader suited for test_detection.py

        Returns a list of simulation steps, reaction names and list of reactant atom lists
        """
        reactions = []
        with open(path) as f:
            lines = f.read().splitlines()
            for line in lines:
                if line[0] == "#" or len(line) == 0:
                    continue
                # we don't care about resids for this
                line, _ = re.subn(r"\([^)]*\)", "", line)
                elems = line.split(";")
                frame, rx = elems[0].split(",")
                # list of atoms
                frags = [
                    [int(atom) for atom in elem.split(",")[2:]] for elem in elems[1:]
                ]
                reactions.append((int(frame), rx, frags))
        return reactions
