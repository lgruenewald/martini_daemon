import warnings
from typing import TextIO

from ..__rust import Fragment
from ..__simulation import Reporter, Simulation


def write_frame(handle: TextIO, sim: Simulation) -> None:
    handle.write(f"==== Frame {sim.current_step} ====\n")
    for frag_id in sim.top.frag_list.get_all_frag_ids():
        frag = sim.top.frag_list.get_fragment(frag_id)
        assert frag is not None
        assert frag.frag_id == frag_id
        atoms: str = ",".join([str(x) for x in frag.atoms])
        handle.write(f"{frag.name},{frag_id};{atoms}\n")
    handle.write("End Frame\n\n")

class FragmentsDump(Reporter):
    def __init__(self) -> None:
        """
        Dump the fragment list, including all atom indices.

        The purpose of this reporter is to facilitate debugging graph matching.

        Uses the file extension ``.fragments_dump``.
        Due to the debug-oriented nature of this format,
        no commitments are made to keep this format forward or backward compatible.
        Due to the debug-oriented nature of this format, it will not be truncated when continuing simulations
        from an older frame.
        """
        warnings.warn("Use the more versatile fragment reporter instead. The use of this class is deprecated.", DeprecationWarning)

    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        self.handle = open(  # noqa: SIM115
            simulation.request_path(".fragments_dump", continue_sim=continue_sim), "a"
        )
        n = simulation.system.num_atoms()
        assert n > 0
        if not continue_sim:
            self.handle.write("# Written by Martini Daemon FragmentsDump\n")
        write_frame(self.handle, simulation)

    def on_simulation_finish(self, simulation: Simulation) -> None:
        self.handle.close()

    def on_reaction(
        self, simulation: Simulation, reactions: list[tuple[str, list[Fragment]]]
    ) -> None:
        write_frame(self.handle, simulation)
