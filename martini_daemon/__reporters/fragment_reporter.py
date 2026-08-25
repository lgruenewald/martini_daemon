import os
import warnings
from collections.abc import Iterator
from dataclasses import dataclass
from typing import TextIO

from ..__rust import Fragment
from ..__simulation import Reporter, Simulation


@dataclass
class FragsFrame:
    step: int
    counts: dict[str, int]
    detailed: bool = False
    frame_index: int | None = None
    time_ps: float | None = None
    fragments: list[Fragment] | None = None

    def serialize(self) -> str:
        """Get the in-file representation of this frame, with a newline at the end."""
        res = []
        res.append(
            f"Step:{self.step},"
            + ",".join([f"{k}:{v}" for k, v in self.counts.items()])
        )
        if self.detailed:
            assert self.frame_index is not None and self.time_ps is not None
            res.append(f"Frame:{self.frame_index},Time:{self.time_ps}")

            assert self.fragments is not None
            for frag in self.fragments:
                frag_name = frag.name
                frag_id = frag.frag_id
                atoms = frag.atoms
                res.append(
                    f"Name:{frag_name},Id:{frag_id},Atoms:[{' '.join(map(str, atoms))}]"
                )

        return "\n".join(res) + "\n"


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
        frags = []
        if self.detailed_frames:
            # skip making this copy if not detailed
            for id in simulation.top.frag_list.get_all_frag_ids():
                frag = simulation.top.frag_list.get_fragment(id)
                # because of type checker
                assert frag is not None
                frags.append(frag)
        frame = FragsFrame(
            simulation.current_step,
            simulation.top.frag_list.frag_counts,
            self.detailed_frames,
            simulation.trajectory_frame,
            simulation.time_ps,
            frags,
        )
        self.handle.write(frame.serialize())

        # self.handle.write(
        #     f"Step:{simulation.current_step},"
        #     + ",".join(
        #         [f"{k}:{v}" for k, v in simulation.top.frag_list.frag_counts.items()]
        #     )
        #     + "\n"
        # )
        # if self.detailed_frames:
        #     self.handle.write(
        #         f"Frame:{simulation.trajectory_frame},Time:{simulation.time_ps}\n"
        #     )
        #     # written in a way to reduce frag_list calls as those go through FFI
        #     for frag_id in simulation.top.frag_list.get_all_frag_ids():
        #         frag = simulation.top.frag_list.get_fragment(frag_id)
        #         assert frag is not None
        #         frag_name = frag.name
        #         frag_atoms = frag.atoms
        #         self.handle.write(
        #             f"Name:{frag_name},Id:{frag_id},Atoms:[{' '.join(map(str, frag_atoms))}]\n"
        #         )

        self.handle.flush()

    def on_simulation_finish(self, simulation: Simulation) -> None:
        assert self.handle is not None
        self.handle.close()

    def interactive_line(self, simulation: Simulation) -> str:
        return f"fragments: {simulation.top.frag_list.num_fragments()}"

    @classmethod
    def iter_fragments(cls, path: str) -> Iterator[FragsFrame]:
        """Read a .frags file frame by frame."""

        class FragsFrameIterator:
            def __init__(self, path: str) -> None:
                self.__handle = open(path)  # noqa: SIM115
                self.__next: str | None = None

            def __advance(self) -> str | None:
                r"""
                Get the next line.

                Returns None on EOF.
                """
                if self.__next is not None:
                    res = self.__next
                    self.__next = None
                    return res
                next = self.__handle.readline()
                if next == "":
                    return None
                return next

            def __backtrack(self, next: str) -> None:
                self.__next = next

            def __finish(self) -> None:
                self.__handle.close()

            def __skip_whitespace(self) -> None:
                """Skip all empty and comment lines."""
                while line := self.__advance():
                    if line is None:
                        break

                    if len(line) > 0 and line[0] != "#":
                        # first content line
                        self.__backtrack(line)
                        break

            def __next__(self) -> FragsFrame:
                self.__skip_whitespace()
                header = self.__advance()
                if header is None:
                    self.__finish()
                    raise StopIteration
                tokens = header.split(",")
                keyword = tokens[0].strip().split(":")[0]
                assert keyword == "Step", ".frags format must start with 'Step:'."
                step = int(tokens[0].strip("Step:"))
                counts = {}
                for tok in tokens[1:]:
                    k, v = tok.strip().split(":")
                    assert k.strip() not in counts
                    counts[k.strip()] = int(v.strip())
                res = FragsFrame(step, counts, False)

                # detailed mode stuff
                while line := self.__advance():
                    if line is None:
                        break
                    line = line.strip()
                    if len(line) == 0 or line[0] == "#":
                        continue
                    tokens = line.split(",")
                    keyword = tokens[0].split(":")[0].strip(":").strip()

                    match keyword:
                        case "Step":
                            # next frame started
                            self.__backtrack(line)
                            break
                        case "Frame":
                            res.detailed = True
                            res.frame_index = int(tokens[0].strip("Frame: "))
                            res.time_ps = float(tokens[1].strip("Time: "))
                        case "Name":
                            res.detailed = True
                            if res.fragments is None:
                                res.fragments = []
                            name = tokens[0].split(":")[1].strip()
                            tok1 = tokens[1].strip()
                            assert tok1.startswith("Id:")
                            id = int(tok1.split(":")[1].strip())
                            tok2 = tokens[2].strip()
                            assert tok2.startswith("Atoms:[") and tok2.endswith("]")
                            atoms = list(map(int, tok2.strip("Atoms:[ ]").split()))
                            res.fragments.append(Fragment(name, id, atoms))
                        case _:
                            warnings.warn(f"Ignoring line '{line}'.")
                return res

            def __iter__(self) -> Iterator[FragsFrame]:
                return self

        return iter(FragsFrameIterator(path))

    @classmethod
    def read_fragments(cls, path: str) -> list[FragsFrame]:
        """
        Read an entire .frags file.

        Note: use iter_fragments if per-frame iteration is desired, as that
        will not load the entire file into memory.
        """
        return list(cls.iter_fragments(path))


class FragCountReporter(FragmentReporter):
    def __init__(self, *args, **kwargs) -> None:  # noqa: ANN002, ANN003
        super().__init__(*args, **kwargs)
        warnings.warn(
            "FragCountReporter has been renamed to FragmentReporter. FragCountReporter is currently aliased, but this is deprecated. ",
            DeprecationWarning,
        )
