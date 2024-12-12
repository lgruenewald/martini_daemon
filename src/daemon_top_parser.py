#!/usr/bin/env python3

"""
A more modular
parser that constructs S* and T* rather than openmm's internal objects
"""

from parser import TopParser, unwrap
from topstar import TopStar, ReactionTemplate
from sysstar import SysStar
import sys
import os
import distutils
import math


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


def DaemonTopFile(file, include_dir=None, defines={}):
    """Parses a Martini Top file for Gromacs and generates T*, sys and top
    from it. Also parses .frag and .rx files included in the .top file.
    """
    # field init
    system = SysStar()
    topology = TopStar(system)

    if include_dir is None:
        include_dir = _get_default_gromacs_include_dir()

    # so .top files stay backwards compatible
    defines["DAEMON"] = 1

    # make parser
    p = TopParser()

    # add all handlers
    def TODO(tokens):
        raise ValueError("Handler not implemented yet")

    def process_defaults(tokens):
        nb_type = unwrap(tokens, 0, "int")
        if nb_type != 1:
            raise ValueError(f"Unsupported non-bonded type: {nb_type}")
        combination_rule = unwrap(tokens, 1, "int")
        if combination_rule != 2:
            raise ValueError("Unsupported combination rule: "
                             f"{combination_rule}")
        if len(tokens) > 2:
            raise ValueError("Too many fields in [ defaults ] directive")

    p.add_level("defaults", process_defaults)

    def process_system(tokens):
        pass
    p.add_level("system", process_system)

    _last_molecule = None

    def process_moltype(tokens):
        nonlocal _last_molecule
        name = unwrap(tokens, 0, "word")
        nrexcl = unwrap(tokens, 1, "int")
        if nrexcl != 1:
            raise ValueError("nrexcl is not 1, only nrexcl=1 is implemented")
        _last_molecule = topology.new_mol_fragment(name)

    p.add_level("moleculetype", process_moltype)

    def process_molecule(tokens):
        name = unwrap(tokens, 0, "word")
        count = unwrap(tokens, 1, "int")
        for i in range(count):
            topology.instantiate(name)

    p.add_level("molecules", process_molecule)

    def last_molecule():
        if last_molecule is None:
            raise ValueError("No [ moleculetype ] given")
        return _last_molecule

    def process_atoms(tokens):
        id = unwrap(tokens, 0, "int") - 1
        type = unwrap(tokens, 1, {"word", "pattern"})
        resnum = unwrap(tokens, 2, "int")
        resname = unwrap(tokens, 3, "word")
        atomname = unwrap(tokens, 4, {"word", "pattern"})
        charge_group_num = unwrap(tokens, 5, "int")
        charge = unwrap(tokens, 6, "float", None)
        mass = unwrap(tokens, 7, "float", None)
        atom_index = len(last_molecule().atoms)
        if id != atom_index:
            raise ValueError("Bad atom ID, are they out of order?"
                             f" got id {id} but expected {atom_index}")
        last_molecule().atoms.append((type, resnum, resname,
                                            atomname, charge_group_num,
                                            charge, mass))

    p.add_level("atoms", process_atoms)

    def process_bonds(tokens):
        i = unwrap(tokens, 0, "int") - 1
        j = unwrap(tokens, 1, "int") - 1
        type = unwrap(tokens, 2, "int")
        if type != 1 and type != 6:
            # 1 is "bond", 6 is "harmonic potential"
            raise ValueError("Unsupported  bond function type")
        length = unwrap(tokens, 3, "float", None)
        force = unwrap(tokens, 4, "float", None)
        last_molecule().bonds.append((system.harmonic_bond, i, j,
                                     [length, force]))
        if type != 6:
            # type 6 does not generate exclusions
            # nrexcl other than 1 is not supported
            last_molecule().add_exclusion(i, j)

    p.add_level("bonds", process_bonds)

    def process_angles(tokens):
        i = unwrap(tokens, 0, "int") - 1
        j = unwrap(tokens, 1, "int") - 1
        k = unwrap(tokens, 2, "int") - 1
        type = unwrap(tokens, 3, "int")
        theta = unwrap(tokens, 4, "float", -1.) * math.pi / 180
        force = unwrap(tokens, 5, "float", -1.)
        match type:
            case 1:
                force_obj = system.harmonic_angle
            case 2:
                force_obj = system.g96_angle
            case 10:
                force_obj = system.restricted_angle
            case _:
                raise ValueError(f"Unsupported angle type {type}.")
        last_molecule().angles.append((force_obj, i, j, k, [theta, force]))

    p.add_level("angles", process_angles)

    def process_dihedrals(tokens):
        i = unwrap(tokens, 0, "int") - 1
        j = unwrap(tokens, 1, "int") - 1
        k = unwrap(tokens, 2, "int") - 1
        l = unwrap(tokens, 3, "int") - 1
        type = unwrap(tokens, 4, "int")
        multiplicity = unwrap(tokens, 7, "int", 1)
        match type:
            case 1 | 9:
                # proper dihedral | proper dihedral (multiple)
                theta = unwrap(tokens, 5, "float")
                force = unwrap(tokens, 6, "float")
                theta_rad = theta * math.pi / 180
                last_molecule().dihedrals.append(
                    (system.proper_dihedral, i, j, k, l,
                     [theta_rad, force, multiplicity])
                )
            case 2:
                # improper
                theta = unwrap(tokens, 5, "float")
                force = unwrap(tokens, 6, "float")
                theta = theta - 360 if theta > 180 else theta
                theta_rad = theta * math.pi / 180
                last_molecule().dihedrals.append(
                    (system.improper_dihedral, i, j, k, l,
                     [theta_rad, force])
                )
            case 3:
                # ryckaert-bellemans = RB
                params = []
                for i in range(6):
                    # C0 to C5
                    params.append(unwrap(tokens, 5+i, "float"))
                last_molecule().dihedrals.append(
                    (system.rb_torsion, i, j, k, l, params)
                )
            case 4:
                # periodic improper dihedral TODO
                raise NotImplementedError
            case 5:
                # Fourier dihedral
                params = []
                for i in range(5):
                    # C1 to C5
                    params.append(unwrap(tokens, 5+i, "float"))
                rb_params = [
                    params[1] + 0.5 * (params[0] + params[2]),
                    0.5 * (-params[0] + 3 * params[2]),
                    -params[1] + 4 * params[3],
                    -2 * params[2],
                    -4 * params[3],
                    0,
                ]
                last_molecule().dihedrals.append(
                    (system.rb_torsion, i, j, k, l, rb_params)
                )
            case 11:
                # combined bending-torsion potential
                force = unwrap(tokens, 5, "float")
                params = [force]
                for i in range(5):
                    # a0 to a4
                    params.append(unwrap(tokens, 6+i, "float"))
                last_molecule().dihedrals.append(
                    (system.combined_bending_torsion, i, j, k, l, params)
                )
            case _:
                raise ValueError("Unsupported dihedral function type {type}.")

    p.add_level("dihedrals", process_dihedrals)

    def process_exclusions(tokens):
        i = unwrap(tokens, 0, "int") - 1
        j = unwrap(tokens, 1, "int") - 1
        last_molecule().add_exclusion(i, j)
        for k in range(2, len(tokens)):
            c = unwrap(tokens, k, "int") - 1
            last_molecule().add_exclusion(i, c)

    p.add_level("exclusions", process_exclusions)

    def process_constraints(tokens):
        i = unwrap(tokens, 0, "int") - 1
        j = unwrap(tokens, 1, "int") - 1
        type = unwrap(tokens, 2, "int")
        length = unwrap(tokens, 3, "float")
        if type == 1 or type == 2:
            last_molecule().constraints.append((i, j, length))
            if type == 1:
                # type 2 doesn't generate exclusions
                last_molecule().add_exclusion(i, j)
        else:
            raise ValueError("Unsupported constraint type")

    p.add_level("constraints", process_constraints)

    p.add_level("pairs", TODO)
    p.add_level("cmap", TODO)

    def process_atomtypes(tokens):
        if len(tokens) != 6:
            raise ValueError("Only atomtypes lines formatted as type, m, q,"
                             "particle_type, V, W are supported.")
        type = unwrap(tokens, 0, "word")
        mass = unwrap(tokens, 1, "float")
        charge = unwrap(tokens, 2, "float")
        particle_type = unwrap(tokens, 3, "word")
        if particle_type != "A":
            raise ValueError("Only A particle type expected in [atomtypes]")
        V = unwrap(tokens, 4, "float")
        W = unwrap(tokens, 5, "float")
        if V != 0.0 or W != 0.0:
            raise ValueError("Only zero V and W are expected in [atomtypes]")
        system.add_atom_type(type, charge, mass)

    p.add_level("atomtypes", process_atomtypes)
    p.add_level("bondtypes", TODO)
    p.add_level("angletypes", TODO)
    p.add_level("dihedraltypes", TODO)
    p.add_level("implicit_genborn_params", TODO)
    p.add_level("pairtypes", TODO)
    p.add_level("cmaptypes", TODO)

    def process_nonbond_params(tokens):
        type1 = unwrap(tokens, 0, "word")
        type2 = unwrap(tokens, 1, "word")
        funct = unwrap(tokens, 2, "int")
        if funct != 1:
            raise ValueError("Only LJ (V/W) non bond params accepted")
        V = unwrap(tokens, 3, "float")
        W = unwrap(tokens, 4, "float")
        system.add_nb_type(type1, type2, V, W)

    p.add_level("nonbond_params", process_nonbond_params)
    p.add_level("virtual_sites2", TODO)
    p.add_level("virtual_sites3", TODO)
    p.add_level("virtual_sitesn", TODO)

    # custom additions: rx and frag
    _last_frag = None

    def last_frag():
        if _last_frag is None:
            raise ValueError("define a [ frag ] first")
        return _last_frag

    def process_frag(tokens):
        nonlocal _last_frag
        name = None
        mol = None
        for i in range(len(tokens)):
            key, value = unwrap(tokens, i, "pair")
            match key:
                case "name":
                    name = value
                case "mol":
                    mol = value
                case _:
                    raise ValueError(f"Unexpected key {key}, "
                                     "expected 'name' or 'mol'")
        if name is None or mol is None:
            raise ValueError("name and mol must be defined in [ frag ]")
        _last_frag = topology.new_frag_fragment(name, mol)

    p.add_level("frag", process_frag)

    def process_fragatoms(tokens):
        parent_id = unwrap(tokens, 0, "int") - 1
        type = unwrap(tokens, 1, {"word", "pattern"})
        name = unwrap(tokens, 2, {"word", "pattern"})
        is_edge = unwrap(tokens, 3, "int") != 0
        last_frag().atoms.append((parent_id, type, name, is_edge))

    p.add_level("fragatoms", process_fragatoms)

    # FIXME this is terrible, the current callback architecture is really
    # unsuitable for this
    last_reaction: ReactionTemplate = None

    def is_complete():
        return (
            last_reaction is not None and last_reaction.r1 is not None and
            last_reaction.r2 is not None and last_reaction.p1 is not None and
            last_reaction.distance_max is not None and
            last_reaction.type is not None
        )

    def process_reaction(tokens):
        nonlocal last_reaction
        for i in range(len(tokens)):
            key, value = unwrap(tokens, i, "pair")
            match key:
                case "name":
                    if last_reaction is not None and not is_complete():
                        raise ValueError(f"unfinished reaction {last_reaction.name}")
                    last_reaction = ReactionTemplate(value, None, None, None, None, None)
                case "r1":
                    last_reaction.r1 = value
                case "r2":
                    last_reaction.r2 = value
                case "p1":
                    last_reaction.p1 = value
                case "type":
                    last_reaction.type = value
                case "distance_max":
                    last_reaction.distance_max = value
            if is_complete():
                topology.new_reaction(last_reaction)
                last_reaction = None

    p.add_level("reaction", process_reaction)

    # run parser
    ok = p.parse(file, include_dir, defines)
    if not ok:
        raise ValueError("Parsing error")

    if last_reaction is not None and not is_complete():
        raise ValueError(f"unfinished reaction {last_reaction.name}")

    return (system, topology)


if __name__ == "__main__":
    argv = sys.argv
    argc = len(argv)
    if argc != 2:
        print("Usage: ./martini_top_parser.py <file>")
        quit(1)
    path = argv[1]

    sys, top = DaemonTopFile(path)
    sys.dump()
    top.dump()
