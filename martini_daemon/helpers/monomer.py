# analysis helper for problems where a group of atoms can be permanently mapped
# into a "monomer" unit, and a graph of monomers

from ..reporters.bond_reporter import read_bonds, _append
import sys
import numpy as np


# convert a .bonds trajectory into a monomer graph based one
def read_monomer_graph(
    path: str, initial_molecules: list[tuple[str, int, int]],
    log_file=sys.stdout, max_bonds=6
):
    """
        Reads a file written by BondReporter. Given a mapping of atoms to
        monomers (initial_molecules), it returns a graph information
        in a similar way as read_bonds does, but for monomers.
        Returns:
        n_frames, n_monomers, max_bonds, empty_value, array
        where array is a dimension of (n_frames, n_monomers, max_bonds) and
        empty_value is used to indicate non values.

        The conversion can take a while, and since the primary
        function of this is to be used in run.py style "recipes"
        this prints progress information to stdout by default. A different
        handle can be used by specifying the log_file argument,
        or can be silenced if it is set to None.

        max_bonds specifies the largest bondedness a monomer can have.

        The result graph can have multiple connections between monomers,
        if there are multiple bonds between the monomers. Bonds, constraints
        and virtual sites are considered bonds, based on how BondReporter
        reports it. In case of vsites, the virtual particle is considered
        bonded to all constructing particles.
    """
    # input data and uninitialized result
    total_monomers = sum([x for _, x, _ in initial_molecules])
    total_frames, total_atoms, bonds_max_bonds, bonds_empty, bonds_data = \
        read_bonds(path)
    monomer_data = np.zeros((total_frames, total_monomers, max_bond), dtype=np.uint32) - 1
    empty = 0xffffffff

    # closure helpers
    def atom_index_to_monomer_index(index: int) -> int | None:
        start = 0
        monomer_index = 0
        for i, (name, n_monomer, n_atoms) in initial_molecules:
            end = start + n_monomer * n_atoms
            if index >= start and index < end:
                return monomer_index + (index - start) // n_atoms
            else:
                start = end
                monomer_index += n_monomer
        return None

    def append(frame, monomer, other_monomer):
        for i in range(max_bonds):
            if monomer_data[frame][monomer][i] == empty:
                monomer_data[frame][monomer][i] = other_monomer
                return
        assert False, f"Too many connections per monomer ({monomer}), max is {max_bonds}."

    # dynamically updated start atom index for every monomer
    for frame in total_frames:
        if frame % 64 == 0 and log_file is not None:
            msg = f"read_monomer_graph: frame {frame} out of {total_frames}."
            if log_file == sys.stdout:
                log_file.write(f"\033[2K\r{msg}")
            else:
                log_file.write(f"{msg}\n")
        current_monomer = 0
        start = 0
        for name, n_monomers, n_atoms in initial_molecules:
            # dynamically updated end atom index for every monomer
            for _ in range(n_monomers):
                end = start + n_atoms
                for atom in range(start, end):
                    # each atom `atom`
                    for other in bonds_data[frame][atom]:
                        # each other atom `other` -> `atom` and `other` are bonded
                        if other == bonds_empty:
                            # no more bonds for atom `atom`
                            break
                        if other >= start and other < end:
                            # `atom` and `other` are the same monomer - uninteresting
                            continue
                        # cross monomer bond, record it
                        other_monomer = atom_index_to_monomer_index(other)
                        append(frame, current_monomer, other_monomer)
                        append(frame, other_monomer, current_monomer)
                start = end
                current_monomer += 1




# calculate degree of conversion
# calculate give #s for different bond orders
# loop analysis - find different types of loops
# chain analysis - find different types of linear chains
