#!/usr/bin/env python3

"""
Refactored https://github.com/maccallumlab/martini_openmm/ into a more modular
parser that constructs S* and T* rather than openmm's internal objects

Also parses fragment types and reaction templates to construct R*
"""

from parser import TopParser, error_at_token, unwrap
from tstar import TopStar, MolFragment
import sys
import os
import distutils


class DaemonTopFile():
    """Parses a Martini Top file for Gromacs and generates T*, sys and top
    from it. Also parses .frag and .rx files included in the .top file.
    """

    def __init__(self, file, include_dir=None, defines={}):
        """
        Load a .top file for martini daemon
        """
        # field init
        self.topology = TopStar([], [], [], [], [])
        self.molecules = []  # convert this to use TopStar()

        # TODO include_dir stuff
        if include_dir is None:
            include_dir = _get_default_gromacs_include_dir()

        # make parser
        p = TopParser()

        # add all handlers
        def TODO(tokens):
            print("Handler not implemented yet")

        p.add_level("defaults", TODO)

        def process_moltype(tokens):
            name = unwrap(tokens, 0, "word")
            nrexcl = unwrap(tokens, 1, "int")
            if nrexcl != 1:
                error_at_token("nrexcl is not 1, only nrexcl=1 is implemented",
                               tokens[1])
            mol_fragment = MolFragment(name, True, [], [], [],
                                       [], [], [], [])
            self.topology.mol_fragments.append(mol_fragment)

        p.add_level("moleculetype", process_moltype)

        def process_molecule(tokens):
            name = unwrap(tokens, 0, "word")
            count = unwrap(tokens, 1, "int")
            self.molecules.append((name, count))

        p.add_level("molecules", process_molecule)

        def last_molecule(tokens):
            if len(self.topology.mol_fragments) == 0:
                error_at_token("No [ moleculetype ] given", tokens[0])
            return self.topology.mol_fragments[-1]

        def process_atoms(tokens):
            id = unwrap(tokens, 0, "int")
            type = unwrap(tokens, 1, "word")
            resnum = unwrap(tokens, 2, "int")
            resname = unwrap(tokens, 3, "word")
            atomname = unwrap(tokens, 4, "word")
            charge_group_num = unwrap(tokens, 5, "int")
            charge = unwrap(tokens, 6, "float", 0.)
            mass = unwrap(tokens, 7, "float", -1.)
            # TODO fix use of special value -1 for "default mass"
            atom_index = len(last_molecule(tokens).atoms) + 1
            if id != atom_index:
                error_at_token("Bad atom ID, are they out of order?"
                               f" got id {id} but expected {atom_index}",
                               tokens[0])
            last_molecule(tokens).atoms.append((type, resnum, resname,
                                                atomname, charge_group_num,
                                                charge, mass))

        p.add_level("atoms", process_atoms)

        def process_bonds(tokens):
            i = unwrap(tokens, 0, "int")
            j = unwrap(tokens, 1, "int")
            type = unwrap(tokens, 2, "type")
            if type != 1 and type != 6:
                # 1 is "bond", 6 is "harmonic potential"
                error_at_token("Unsupported  bond function type", tokens[0])
            length = unwrap(tokens, 3, "float", -1.)
            force = unwrap(tokens, 4, "float", -1.)
            last_molecule(tokens).harmonic_bonds.append((i, j, length, force))
            if type != 6:
                # type 6 does not generate exclusions
                # nrexcl other than 1 is not supported
                last_molecule(tokens).add_exclusion(i, j)

        p.add_level("bonds", process_bonds)

        def process_angles(tokens):
            i = unwrap(tokens, 0, "int")
            j = unwrap(tokens, 1, "int")
            k = unwrap(tokens, 2, "int")
            type = unwrap(tokens, 3, "int")
            theta = unwrap(tokens, 4, "float", -1.)
            force = unwrap(tokens, 5, "float", -1.)
            if type == 1:
                # harmonic
                last_molecule(tokens).harmonic_angles.append((i, j, k,
                                                             theta, force))
            else:
                # TODO 2 is g96 angle
                # 10 is restricted angle, they are also supported by
                # martini_openmm
                error_at_token("Unsupported angle function type", tokens[0])

        p.add_level("angles", process_angles)

        def process_dihedrals(tokens):
            i = unwrap(tokens, 0, "int")
            j = unwrap(tokens, 1, "int")
            k = unwrap(tokens, 2, "int")
            l = unwrap(tokens, 3, "int")
            type = unwrap(tokens, 4, "int")
            theta = unwrap(tokens, 5, "float", -1.)
            force = unwrap(tokens, 6, "float", -1.)
            multiplicity = unwrap(tokens, 7, "int", 1)
            if type == 1:
                # proper dihedral
                last_molecule(tokens).proper_dihedrals.append(
                    (i, j, k, l, theta, force, multiplicity)
                )
            elif type == 2:
                # improper
                last_molecule(tokens).improper_dihedrals.append(
                    (i, j, k, l, theta, force)
                )
            else:
                # TODO 2, 3, 4, 5, 9 and 11 should also supported by
                # martini_openmm
                error_at_token("Unsupported dihedral function type", tokens[0])

        p.add_level("dihedrals", process_dihedrals)

        def process_exclusions(tokens):
            i = unwrap(tokens, 0, "int")
            j = unwrap(tokens, 1, "int")
            last_molecule(tokens).add_exclusion(i, j)
            for i in range(2, len(tokens)):
                c = unwrap(tokens, i, "int")
                last_molecule(tokens).add_exclusion(i, c)

        p.add_level("exclusions", process_exclusions)

        def process_constraints(tokens):
            i = unwrap(tokens, 0, "int")
            j = unwrap(tokens, 1, "int")
            type = unwrap(tokens, 2, "int")
            length = unwrap(tokens, 3, "float")
            if type == 1 or type == 2:
                last_molecule(tokens).constraints.append((i, j, length))
                if type == 1:
                    # type 2 doesn't generate exclusions
                    last_molecule(tokens).add_exclusion(i, j)
            else:
                error_at_token("Unsupported constraint type", tokens[0])

        p.add_level("constraints", process_constraints)

        p.add_level("pairs", TODO)
        p.add_level("cmap", TODO)
        p.add_level("atomtypes", TODO)
        p.add_level("bondtypes", TODO)
        p.add_level("angletypes", TODO)
        p.add_level("dihedraltypes", TODO)
        p.add_level("implicit_genborn_params", TODO)
        p.add_level("pairtypes", TODO)
        p.add_level("cmaptypes", TODO)
        p.add_level("nonbond_params", TODO)
        p.add_level("virtual_sites2", TODO)
        p.add_level("virtual_sites3", TODO)
        p.add_level("virtual_sitesn", TODO)

        # run parser
        p.parse(file, include_dir, defines)


if __name__ == "__main__":
    argv = sys.argv
    argc = len(argv)
    if argc != 2:
        print("Usage: ./martini_top_parser.py <file>")
        quit(1)
    path = argv[1]

    daemontop = DaemonTopFile(path)
    print("DUMP OF CREATED Topology, T* and System")
    print("==== Topology / Molecule Fragment Types ====")
    for molfrag in daemontop.topology.mol_fragments:
        print(molfrag)


def _get_default_gromacs_include_dir():
    """Find the location where gromacs #include files are referenced from, by
    searching for (1) gromacs environment variables, (2) for the gromacs binary
    'pdb2gmx' or 'gmx' in the PATH, or (3) just using the default gromacs
    install location, /usr/local/gromacs/share/gromacs/top

    Directly taken from openmm GromacsTopParser
    https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/gromacstopfile.py
    """
    if "GMXDATA" in os.environ:
        return os.path.join(os.environ["GMXDATA"], "top")
    if "GMXBIN" in os.environ:
        return os.path.abspath(
            os.path.join(os.environ["GMXBIN"], "..", "share", "gromacs", "top")
        )

    pdb2gmx_path = distutils.spawn.find_executable("pdb2gmx")
    if pdb2gmx_path is not None:
        return os.path.abspath(
            os.path.join(os.path.dirname(pdb2gmx_path), "..", "share",
                         "gromacs", "top")
        )
    else:
        gmx_path = distutils.spawn.find_executable("gmx")
        if gmx_path is not None:
            return os.path.abspath(
                os.path.join(os.path.dirname(gmx_path), "..", "share",
                             "gromacs", "top")
            )

    return "/usr/local/gromacs/share/gromacs/top"
