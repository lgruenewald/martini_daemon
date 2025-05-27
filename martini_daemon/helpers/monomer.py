# analysis helper for problems where a group of atoms can be permanently mapped
# into a "monomer" unit, and a graph of monomers

from ..reporters.bond_reporter import read_bonds
import sys
import numpy as np


# convert a .bonds trajectory into a monomer graph based one
def read_monomer_graph(
    path: str, initial_molecules: list[tuple[str, int, int]]
):
    """
        Reads a file written by BondReporter. Given a mapping of atoms to
        monomers (initial_molecules), it returns a graph information
        in a similar way as read_bonds does, but for monomers.
        Returns: frames - a list of frames, each frame a numpy array
        of dimension (n_bonds, 2), but in this case it is between
        monomers, not between atoms.

        The conversion can take a while, and since the primary
        function of this is to be used in run.py style "recipes"
        this prints progress information to stdout by default.

        The result graph can have multiple connections between monomers,
        if there are multiple bonds between the monomers. Bonds, constraints
        and virtual sites are considered bonds, based on how BondReporter
        reports it. In case of vsites, the virtual particle is considered
        bonded to all constructing particles.
    """
    # input data and uninitialized result
    # make a mapping atom -> monomer
    total_atoms = sum([x * y for _, x, y in initial_molecules])
    mapping = np.empty(total_atoms, np.uint32)
    atom = 0
    monomer = 0
    for _, monomer_count, atom_count in initial_molecules:
        for _ in range(monomer_count):
            mapping[atom:atom+atom_count] = monomer
            atom += atom_count
            monomer += 1
    bond_frames = read_bonds(path)
    monomer_frames = []
    for frame in bond_frames:
        # bonds that are between different monomers
        cross_bonds = filter(
            lambda pair: pair[0] != pair[1],
            map(
                lambda pair: (mapping[pair[0]], mapping[pair[1]]),
                frame
            )
        )
        # create a np array for them
        monomer_frames.append(np.fromiter(
            cross_bonds,
            dtype=np.dtype((np.uint32, 2))
        ))
    return monomer_frames




# calculate degree of conversion
# calculate give #s for different bond orders
# loop analysis - find different types of loops
# chain analysis - find different types of linear chains
