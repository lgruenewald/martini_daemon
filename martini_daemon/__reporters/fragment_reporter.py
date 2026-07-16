import os
import warnings
from dataclasses import dataclass
from typing import TextIO

from ..__rust import Fragment, PeriodicBox
from ..__simulation import Reporter, Simulation


@dataclass
class FragsFrame:
    step: int
    counts: dict[str, int]
    frame_index: int | None = None
    time_ps: float | None = None
    box: PeriodicBox | None = None
    fragments: list[Fragment] | None = None


class FragmentReporter(Reporter):
    def __init__(
        self,
        detailed: bool = False,
        on_trajectory: bool = True,
        on_reaction: bool = False,
    ) -> None:
        """
        Write fragment information per trajectory frame to a `.frags` file.

        Without the detailed flag, only simulation steps and fragment counts are written, useful for rate checks.

        With the detailed_frames flag, additional frame metadata is written, alongside a full list of frags.
        Useful for querying, but takes more disk space.

        With the on_trajectory flag, fragment info is written on trajectory frames. By default, this is on.

        With the on_reaction flag, fragment info is written when reactions occur.
        This can reduce disk space usage in some scenarios if on_trajectory is off,
        since only reactions change fragment counts, but might make analysis scripts more convoluted to write.
        If there are more reactions than trajectory frames, or if on_trajectory is not turned off,
        this will increase disk space usage.
        """
        self.detailed_frames = detailed
        self.write_on_trajectory = on_trajectory
        self.write_on_reaction = on_reaction
        self.path: str | None = None
        self.handle: TextIO | None = None

    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        self.path = simulation.request_path(".frags", continue_sim=continue_sim)
        truncate = continue_sim and os.path.exists(self.path)
        self.handle = open(self.path, "r+" if truncate else "w")  # noqa: SIM115

        if truncate:
            truncate_to = simulation.current_step
            tell = self.handle.tell()
            while (line := self.handle.readline()) != "":
                tokens = line.split(",")
                token0 = tokens[0].strip()
                if token0.startswith("Step"):
                    c_step = int(token0.strip("Step:"))
                    if c_step > truncate_to:
                        break
                tell = self.handle.tell()
            self.handle.truncate(tell)
            self.handle.seek(0, os.SEEK_END)

    def on_trajectory_frame(self, simulation: Simulation) -> None:
        if self.write_on_trajectory:
            self.__write_frame(simulation)

    def on_reaction(
        self, simulation: Simulation, reactions: list[tuple[str, list[Fragment]]]
    ) -> None:
        if self.write_on_reaction:
            self.__write_frame(simulation)

    def __write_frame(self, simulation: Simulation) -> None:
        assert self.handle is not None
        self.handle.write(
            f"Step:{simulation.current_step},"
            + ",".join(
                [f"{k}:{v}" for k, v in simulation.top.frag_list.frag_counts.items()]
            )
            + "\n"
        )
        if self.detailed_frames:
            self.handle.write(
                f"Frame:{simulation.trajectory_frame},Time:{simulation.time_ps}\n"
            )
            # written in a way to reduce frag_list calls as those go through FFI
            for frag_id in simulation.top.frag_list.get_all_frag_ids():
                frag = simulation.top.frag_list.get_fragment(frag_id)
                assert frag is not None
                frag_name = frag.name
                frag_atoms = frag.atoms
                self.handle.write(
                    f"Name:{frag_name},Id:{frag_id},Atoms:[{' '.join(map(str, frag_atoms))}]\n"
                )

        self.handle.flush()

    def on_simulation_finish(self, simulation: Simulation) -> None:
        assert self.handle is not None
        self.handle.close()

    def interactive_line(self, simulation: Simulation) -> str:
        return f"fragments: {simulation.top.frag_list.num_fragments()}"

    @classmethod
    def read_fragments(cls, path: str) -> list[FragsFrame]:
        res: list[FragsFrame] = []
        with open(path) as handle:
            lines = handle.readlines()

        for line in lines:
            tokens = line.split(",")
            keyword = tokens[0].strip().split(":")[0]
            match keyword:
                case "Step":
                    # frame header
                    step = int(tokens[0].strip("Step:"))
                    counts = {}
                    for tok in tokens[1:]:
                        k, v = tok.strip().split(":")
                        assert k.strip() not in counts
                        counts[k.strip()] = int(v.strip())
                    res.append(FragsFrame(step, counts))
                case "Frame":
                    # secondary frame header if detailed mode is on
                    assert len(res) > 0, (
                        ".frags file corrupt, as new frames must start with a line Step:... first."
                    )
                    res[-1].frame_index = int(tokens[0].strip("Frame: "))
                    res[-1].time_ps = float(tokens[1].strip("Time: "))
                    # if the Frame line is present => detailed mode was on when writing this
                    res[-1].fragments = []
                case "Name":
                    # a single fragment entry, only written if detailed mode is on
                    if res[-1].fragments is None:
                        # normally shouldn't occur on files written by this reporter
                        res[-1].fragments = []
                    name = tokens[0].split(":")[1].strip()
                    tok1 = tokens[1].strip()
                    assert tok1.startswith("Id:")
                    id = int(tok1.split(":")[1].strip())
                    tok2 = tokens[2].strip()
                    assert tok2.startswith("Atoms:[") and tok2.endswith("]")
                    atoms = list(map(int, tok2.strip("Atoms:[ ]").split()))
                    res[-1].fragments.append(Fragment(name, id, atoms))
                # forward compat (?): ignore lines not starting with a known keyword

        return res


class FragCountReporter(FragmentReporter):
    def __init__(self, *args, **kwargs) -> None:  # noqa: ANN002, ANN003
        super().__init__(*args, **kwargs)
        warnings.warn(
            "FragCountReporter has been renamed to FragmentReporter. FragCountReporter is currently aliased, but this is deprecated. ",
            DeprecationWarning,
        )
