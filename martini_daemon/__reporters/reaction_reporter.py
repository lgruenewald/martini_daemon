import os
import re
from typing import NamedTuple, TextIO

from ..__formats import TrajectoryWriter
from ..__rust import Fragment
from ..__simulation import Reporter, Simulation


def truncate_reactions(handle: TextIO, current_step: int) -> None:
    prev_pos = 0
    while len(line := handle.readline()) > 0:
        if len(line) > 0 and line[0] != "#":
            frame = int(line.split(",")[0])
            if frame > current_step:
                break
        prev_pos = handle.tell()
    handle.truncate(prev_pos)
    handle.seek(0, os.SEEK_END)

class ReactionsFileReactant(NamedTuple):
    """Reactant parsed from a .reactions file."""
    # FIXME: identical signature to Fragment, so this is redundant, but for backwards compat it's a named tuple
    frag_name: str
    frag_id: int
    atoms: list[int]

class ReactionsFileReaction(NamedTuple):
    """Reaction parsed from a .reactions file."""
    # FIXME: normally I would prefer a dataclass, but a named tuple is used for backwards compat
    step: int
    reaction_name: str
    reactant: list[ReactionsFileReactant]

class ReactionReporter(Reporter):
    def __init__(self, traj_format: str | None = None) -> None:
        """
        Create a ReactionReporter to log reaction information to <name>.reactions.

        The created file has a text format, where every line is a reaction.
        First, the frame number and reaction name are separated by a comma,
        then, each reactant is separated by a semicolon. Each reactant
        will have its frag name, internal frag id, and atom indices
        (-1 for missing optional or forbidden atoms) printed. Atom indices are 0 indexed.

        The residue IDs involved are also included as a comma separated list, within parentheses, separately
        for each reactant.

        Example:
        frame,reaction_name;reactant1_name,reactant1_id(res:resid1,...residn),atoms...;...reactantn_name(res:resid1,...residn),reactantn_id,atoms...

        If traj_format is set to a string, a trajectory file of that format will be created, writing a frame for all MD
        steps where reactions occured (before minimization).

        """
        # .reactions
        self.path: str | None = None
        self.handle: TextIO | None = None
        # _reactions.{format}
        self.traj_format = traj_format if traj_format is None or traj_format[0] == "." else f".{traj_format}"
        self.traj_path: str | None = None
        self.traj_writer: TrajectoryWriter | None = None

    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        # .reactions
        self.path = simulation.request_path(".reactions", continue_sim=continue_sim)
        truncate = continue_sim and os.path.exists(self.path)
        self.handle = open(self.path, "r+" if truncate else "w")  # noqa: SIM115

        if truncate:
            truncate_reactions(self.handle, simulation.current_step)
        else:
            self.handle.write(
                "# sim step,reaction_name;"
                + "reactant1_name,reactant1_id(res:resid1,...residn),atoms...;..."
                + "reactantn_name,reactantn_id(res:resid1,...residn),atoms...;\n",
            )

        # _reactions.{format}
        if self.traj_format is not None:
            self.traj_path = simulation.request_path(f"_reactions{self.traj_format}", continue_sim=continue_sim)
            append = continue_sim and os.path.exists(self.traj_path)
            self.traj_writer = TrajectoryWriter(
                self.traj_path,
                append=append,
                # formats have varying metadata on sim time / step, so truncate based on the # of frames written before
                # sim.trajectory_frame is the number of frames that were finished writing
                # (or during trajectory frame also the index of the frame currently being written)
                keep_n_frames=simulation.trajectory_frame if append else None,
            )


    def on_simulation_finish(self, simulation: Simulation) -> None:
        assert self.handle is not None
        self.handle.close()
        self.handle = None
        if self.traj_writer is not None:
            self.traj_writer.close()
            self.traj_writer = None

    @classmethod
    def __get_resids(cls, sim: Simulation, atoms: list[int]) -> str:
        res = set()
        for atom in atoms:
            if atom != -1:
                res.add(f"{sim.system.get_res_ids()[atom]}")
        return "(res:" + ",".join(res) + ")"

    # need to be post, as the modification algorithm can reject some reactions
    def on_reaction(
        self, simulation: Simulation, reactions: list[tuple[str, list[Fragment]]]
    ) -> None:
        # .reactions
        assert self.handle is not None
        for rx, frags in reactions:
            self.handle.write(
                f"{simulation.current_step},{rx};"
                + ";".join(
                    [
                        f"{frag.name},{frag.frag_id}"
                        f"{self.__get_resids(simulation, frag.atoms)},"
                        + ",".join([f"{atom}" for atom in frag.atoms])
                        for frag in frags
                    ]
                )
                + "\n",
            )
        self.handle.flush()
        # _reactions.{format}
        if self.traj_format is not None:
            assert self.traj_writer is not None
            pos, box = simulation.context.get_positions()
            self.traj_writer.write_frame(simulation.current_step, simulation.time_ps, box, pos)

    def interactive_line(self, simulation: Simulation) -> str:
        return f"reactions: {simulation.reactions_so_far}"

    @classmethod
    def read_reactions(
        cls, path: str
    ) -> list[ReactionsFileReaction]:
        """
        .reactions format reader suited for test_detection.py.

        Returns a list of simulation steps, reaction names and list of reactant atom lists
        """
        reactions: list[ReactionsFileReaction] = []
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
                frags: list[ReactionsFileReactant] = []
                for elem in elems[1:]:
                    tokens = elem.split(",")
                    # name, frag_id, atoms
                    frags.append(
                        ReactionsFileReactant(tokens[0], int(tokens[1]), [int(tok) for tok in tokens[2:]])
                    )
                reactions.append(ReactionsFileReaction(int(frame), rx, frags))
        return reactions
