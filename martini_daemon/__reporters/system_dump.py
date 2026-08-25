import re
from typing import TextIO

from ..__core import BondedForce
from ..__rust import Fragment
from ..__simulation import Reporter, Simulation


def write_frame(handle: TextIO, sim: Simulation) -> None:
    handle.write(f"==== Frame {sim.current_step} ====\n")

    sys = sim.system
    handle.write("Atoms\n")
    handle.write("# (name, resid, resname, type, charge, mass, sc_lam, sc_alpha)\n")

    handle.writelines(
        f"('{sys.get_name(atom)}', {sys.get_res_id(atom)}, '{sys.get_res_name(atom)}', "
        + f"'{sys.get_type(atom)}', {sys.get_charge(atom)}, {sys.get_mass(atom)}, "
        + f"{sys.get_sc_lam(atom)}, {sys.get_sc_alpha(atom)})\n"
        for atom in range(sys.num_atoms())
    )

    handle.write("\n")
    handle.write("Forces\n")
    for force in sys.get_forces():
        if not isinstance(force, BondedForce):
            continue
        if len([force.iterate_bonds()]) == 0:
            continue
        handle.write(f"Force:{force.get_name()}\n")
        handle.writelines(
            "(" + ", ".join(str(x) for x in members + params) + ")\n"
            for _, (members, params) in force.iterate_bonds()
        )

    handle.write("End Frame\n\n")


class SystemDump(Reporter):
    def __init__(self) -> None:
        """
        Create a SystemDump Reporter.

        Dumps all info from System, including all atom details and all interactions
        in a human-readable plaintext file. Dumps it at the start of a simulation and
        when there is any reactions.

        The purpose of this reporter is twofold:

        * to facilitate debugging reactions,
        * it is also used in the modification algorithm unit tests.

        Uses the file extension ``.system_dump``.
        Due to the debug-oriented nature of this format,
        no commitments are made to keep this format forward or backward compatible.
        Due to the debug-oriented nature of this format, it will not be truncated when continuing simulations
        from an older frame.
        """

    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        self.handle = open(  # noqa: SIM115
            simulation.request_path(".system_dump", continue_sim=continue_sim), "a"
        )
        n = simulation.system.num_atoms()
        assert n > 0
        if not continue_sim:
            self.handle.write("# Written by Martini Daemon SystemDump\n")
        write_frame(self.handle, simulation)

    def on_simulation_finish(self, simulation: Simulation) -> None:
        self.handle.close()

    def on_reaction(
        self, simulation: Simulation, reactions: list[tuple[str, list[Fragment]]]
    ) -> None:
        write_frame(self.handle, simulation)

    @classmethod
    def read_dump(
        cls,
        path: str,
    ) -> list[
        tuple[
            int,
            list[tuple[str, int, str, str, float, float, float, float]],
            list[tuple[str, list[list[float]]]],
        ]
    ]:
        """
        .sstar dump reader.

        Returns a list of frames read.

        * Each frame is a tuple of sim step, atoms and forces.
        * atoms is a list of tuples of:
            * name
            * resid
            * resname
            * atom type
            * charge
            * mass
            * two additional float parameters for soft core
        * forces is a list of tuples of:
            * force name
            * interaction list, which is a list of tuples.
                * these tuples contain members and parameters, with count specific to the interaction at hand (floats).
        """
        frames = []
        with open(path) as f:
            lines = f.read().splitlines()

        prev_i = None
        i = 0

        def advance() -> str | None:
            # return the next non-comment, non-empty line in file, or None if we reached the end of file
            nonlocal i, prev_i
            prev_i = i
            while i < len(lines):
                i += 1
                if len(lines[i - 1]) > 0 and lines[i - 1][0] != "#":
                    return lines[i - 1]
            return None

        def backtrack() -> None:
            nonlocal i, prev_i
            assert prev_i is not None
            i = prev_i
            prev_i = None

        def parse_atoms(
            atoms: list[tuple[str, int, str, str, float, float, float, float]],
        ) -> None:

            while line := advance():
                if line[0] != "(":
                    backtrack()
                    break
                name, res_id, res_name, atom_type, charge, mass, sc_lam, sc_alpha = (
                    line.strip("()").replace(" ", "").split(",")
                )
                atoms.append(
                    (
                        name,
                        int(res_id),
                        res_name,
                        atom_type,
                        float(charge),
                        float(mass),
                        float(sc_lam),
                        float(sc_alpha),
                    )
                )

        def parse_force(force: list[list[float]]) -> None:
            while line := advance():
                if line == "None":
                    continue
                if line[0] != "(":
                    backtrack()
                    break
                elems = [float(x) for x in line.strip("()").replace(" ", "").split(",")]
                force.append(elems)

        def parse_forces(forces: list[tuple[str, list[list[float]]]]) -> None:
            while line := advance():
                if line[:6] != "Force:":
                    backtrack()
                    break
                force = []
                forces.append((line[6:], force))
                parse_force(force)

        def parse_frame(line: str) -> None:
            nonlocal frames
            pat = re.compile("[0-9]+")
            num = pat.search(line)
            assert num is not None
            atoms = []
            forces = []
            frame = (int(line[num.start() : num.end()]), atoms, forces)
            frames.append(frame)

            while cline := advance():
                match cline:
                    case "Atoms":
                        parse_atoms(atoms)
                    case "Forces":
                        parse_forces(forces)
                    case "End Frame":
                        break
                    case _:
                        assert False, f"Unexpected line: {cline}"

        while cline := advance():
            if re.match(r"^==== Frame [0-9]+ ====$", cline):
                parse_frame(cline)

        return frames
