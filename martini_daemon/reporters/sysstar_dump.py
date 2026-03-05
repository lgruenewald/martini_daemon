from .reporter import Reporter
import re

class SysStarDump(Reporter):
    """
        Dumps all info from S* -- all atom details and all interactions
        in a human-readable plaintext file. Dumps it at the start and
        when there is any change.
    """

    def __init__(self):
        pass

    def on_set_xtc_path(self, xtc_name):
        self._open(xtc_name + ".sstar")
        n = self._sysstar.len_atoms()
        assert n > 0
        self._print("Format: SStar Dump")
        self._print("Version: 0")
        self._print(f"N: {n}")
        self._print("# Written by Martini Daemon SysStarDump")
        self._write_frame(0)

    def post_modification(self, i, name, reactions):
        self._write_frame(i)

    def _write_frame(self, i):
        self._print(f"==== Frame {i} ====")

        sys = self._sysstar
        self._print("Atoms")
        self._print("# (name, resid, resname, type, charge, mass, sc_lam, sc_alpha)")
        for deets in sys._atom_list:
            self._print(f"{deets}")

        self._print("")
        self._print("Forces")
        for force in sys.modular_forces:
            if force._list is None or len(force._list) == 0:
                continue
            else:
                self._print(f"Force:{type(force).__name__}")
                for deets in force._list:
                    self._print(f"{deets}")

        self._print("End Frame")
        self._print("")

    @staticmethod
    def read_dump(path: str):
        """
        .sstar dump reader
        TODO the format writer and reader should go together in the same place
        """
        frames = []
        with open(path, "r") as f:
            lines = f.read().splitlines()

        prev_i = None
        i = 0

        def advance():
            # return the next non-comment, non-empty line in file, or None if we reached the end of file
            nonlocal i, prev_i
            prev_i = i
            while i < len(lines):
                i += 1
                if len(lines[i-1]) > 0 and lines[i-1][0] != "#":
                    return lines[i-1]
            return None

        def backtrack():
            nonlocal i, prev_i
            assert prev_i is not None
            i = prev_i
            prev_i = None


        # TODO header parsing

        def parse_atoms(atoms):
            while line := advance():
                if line[0] != "(":
                    backtrack()
                    break
                name, res_id, res_name, atom_type, charge, mass, *other = line.strip("()").replace(" ", "").split(",")
                atoms.append(
                    (name, int(res_id), res_name, atom_type, float(charge), float(mass), *other)
                )

        def parse_force(force):
            while line := advance():
                if line == "None":
                    continue
                if line[0] != "(":
                    backtrack()
                    break
                elems = [
                    float(x) for x in line.strip("()").replace(" ", "").split(",")
                ]
                force.append(elems)

        def parse_forces(forces):
            while line := advance():
                if line[:6] != "Force:":
                    backtrack()
                    break
                force = []
                forces.append((line[6:], force))
                parse_force(force)

        def parse_frame(line):
            nonlocal frames
            pat = re.compile("[0-9]+")
            num = pat.search(line)
            atoms = []
            forces = []
            frame = (int(line[num.start():num.end()]), atoms, forces)
            frames.append(frame)

            while line := advance():
                match line:
                    case "Atoms":
                        parse_atoms(atoms)
                    case "Forces":
                        parse_forces(forces)
                    case "End Frame":
                        break
                    case _:
                        assert False, f"Unexpected line: {line}"

        while cline := advance():
            if re.match(r"^==== Frame [0-9]+ ====$", cline):
                parse_frame(cline)

        return frames