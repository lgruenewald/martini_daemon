#!/usr/bin/env python3

"""
A more modular
parser that constructs S* and T* rather than openmm's internal objects
"""

from .parser import TopParser, unwrap
from .topstar import TopStar, ReactionTemplate
from .sysstar import SysStar
import sys
import os
import distutils
import math
from openmm.unit import nanometer


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


def DaemonTopFile(file, include_dir=None, defines={},
                  epsilon_r=15.0, nonbonded_cutoff=1.1*nanometer):
    """Parses a Martini Top file for Gromacs and generates T*, sys and top
    from it. Also parses .frag and .rx files included in the .top file.
    """
    # field init
    system = SysStar(epsilon_r, nonbonded_cutoff)
    topology = TopStar(system)

    if include_dir is None:
        include_dir = _get_default_gromacs_include_dir()

    # so .top files stay backwards compatible
    defines["DAEMON"] = 1

    # make parser
    p = TopParser()

    system_defined = False  # whether the [system] happened yet

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
        id = unwrap(tokens, 0, "index")
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
        i = unwrap(tokens, 0, "index")
        j = unwrap(tokens, 1, "index")
        type = unwrap(tokens, 2, "int")
        length = unwrap(tokens, 3, "float")
        if type == 1 or type == 6:
            # harmonic bond / harmonic potential
            force = unwrap(tokens, 4, "float")
            last_molecule().bonds.append((system.harmonic_bond, i, j,
                                         [length, force]))
        elif type == 3:
            # morse
            D = unwrap(tokens, 4, "float")
            beta = unwrap(tokens, 5, "float")
            last_molecule().bonds.append((system.morse_bond, i, j,
                                         [length, D, beta]))
        elif type == 4:
            # cubic bond
            kb = unwrap(tokens, 4, "float")
            kcub = unwrap(tokens, 5, "float")
            last_molecule().bonds.append((system.cubic_bond, i, j,
                                         [length, kb, kcub]))
        elif type == 7:
            # FENE (finitely extensible nonlinear elastic) bond
            force = unwrap(tokens, 4, "float")
            last_molecule().bonds.append((system.fene_bond, i, j,
                                         [length, force]))
        else:
            raise ValueError("Unsupported  bond function type")
        if type != 6:
            # type 6 does not generate exclusions
            # nrexcl other than 1 is not supported
            last_molecule().add_exclusion(i, j)

    p.add_level("bonds", process_bonds)

    def process_angles(tokens):
        i = unwrap(tokens, 0, "index")
        j = unwrap(tokens, 1, "index")
        k = unwrap(tokens, 2, "index")
        type = unwrap(tokens, 3, "int")
        theta = unwrap(tokens, 4, "float") * math.pi / 180
        force = unwrap(tokens, 5, "float")
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
        i = unwrap(tokens, 0, "index")
        j = unwrap(tokens, 1, "index")
        k = unwrap(tokens, 2, "index")
        l = unwrap(tokens, 3, "index")
        type = unwrap(tokens, 4, "int")
        match type:
            case 1 | 9:
                # proper dihedral | proper dihedral (multiple)
                theta = unwrap(tokens, 5, "float")
                force = unwrap(tokens, 6, "float")
                theta_rad = theta * math.pi / 180
                multiplicity = unwrap(tokens, 7, "int", 1)
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
                for c in range(6):
                    # C0 to C5
                    params.append(unwrap(tokens, 5+c, "float"))
                last_molecule().dihedrals.append(
                    (system.rb_torsion, i, j, k, l, params)
                )
            case 4:
                # periodic improper dihedral TODO
                raise NotImplementedError
            case 5:
                # Fourier dihedral
                params = []
                for c in range(4):
                    # C1 to C4
                    params.append(unwrap(tokens, 5+c, "float"))
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
                for c in range(5):
                    # a0 to a4
                    params.append(unwrap(tokens, 6+c, "float"))
                last_molecule().dihedrals.append(
                    (system.combined_bending_torsion, i, j, k, l, params)
                )
            case _:
                raise ValueError("Unsupported dihedral function type {type}.")

    p.add_level("dihedrals", process_dihedrals)

    def process_exclusions(tokens):
        i = unwrap(tokens, 0, "index")
        j = unwrap(tokens, 1, "index")
        last_molecule().add_exclusion(i, j)
        for k in range(2, len(tokens)):
            c = unwrap(tokens, k, "index")
            last_molecule().add_exclusion(i, c)

    p.add_level("exclusions", process_exclusions)

    def process_constraints(tokens):
        i = unwrap(tokens, 0, "index")
        j = unwrap(tokens, 1, "index")
        type = unwrap(tokens, 2, "int")
        length = unwrap(tokens, 3, "float")
        if type == 1 or type == 2:
            last_molecule().bonds.append((system.constraint, i, j, [length]))
            if type == 1:
                # type 2 doesn't generate exclusions
                last_molecule().add_exclusion(i, j)
        else:
            raise ValueError("Unsupported constraint type")

    p.add_level("constraints", process_constraints)

    def process_pairtypes(tokens):
        t1 = unwrap(tokens, 0, "word")
        t2 = unwrap(tokens, 0, "word")
        type = unwrap(tokens, 2, "int")
        if type != 1:
            raise ValueError("Unsupported pairs type")
        V = unwrap(tokens, 3, "float")
        W = unwrap(tokens, 4, "float")
        system.pairs.add_type(t1, t2, V, W)

    p.add_level("pairtypes", process_pairtypes)

    def process_pairs(tokens):
        i = unwrap(tokens, 0, "index")
        j = unwrap(tokens, 1, "index")
        type = unwrap(tokens, 2, "int")
        if type != 1:
            raise ValueError("Unsupported pairs type")
        if len(tokens) >= 5:
            V = unwrap(tokens, 3, "float")
            W = unwrap(tokens, 4, "float")
            last_molecule().interactions.append(
                (system.pairs, [i, j], [V, W])
            )
        else:
            last_molecule().interactions.append((system.pairs, [i, j], []))

    p.add_level("pairs", process_pairs)
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

    def process_virtual_sites1(tokens):
        vid = unwrap(tokens, 0, "index")
        member = unwrap(tokens, 1, "index")
        type = unwrap(tokens, 2, "int")
        if type != 1:
            raise ValueError(f"Virtual site 1 type {type} not implemented.")
        last_molecule().interactions.append(
            (system.vsite_avg, [vid, member], [1.])
        )

    p.add_level("virtual_sites1", process_virtual_sites1)

    def process_virtual_sites2(tokens):
        vid = unwrap(tokens, 0, "index")
        i = unwrap(tokens, 1, "index")
        j = unwrap(tokens, 2, "index")
        members = [vid, i, j]
        type = unwrap(tokens, 3, "int")
        match type:
            case 1:
                a = unwrap(tokens, 4, "float")
                weights = [1-a, a]
                last_molecule().interactions.append(
                    (system.vsite_avg, members, weights)
                )
            case 2:
                d = unwrap(tokens, 4, "float")
                last_molecule().interactions.append(
                    (system.vsite_2fd, members, [d])
                )
            case _:
                raise ValueError(
                    f"Virtua site 2 type {type} not implemented."
                )

    p.add_level("virtual_sites2", process_virtual_sites2)

    def process_virtual_sites3(tokens):
        vid = unwrap(tokens, 0, "index")
        i = unwrap(tokens, 1, "index")
        j = unwrap(tokens, 2, "index")
        k = unwrap(tokens, 3, "index")
        members = [vid, i, j, k]
        type = unwrap(tokens, 4, "int")
        match type:
            case 1:
                # 3
                a = unwrap(tokens, 5, "float")
                b = unwrap(tokens, 6, "float")
                weights = [1 - a - b, a, b]
                last_molecule().interactions.append(
                    (system.vsite_avg, members, weights)
                )
            case 2:
                # 3fd
                a = unwrap(tokens, 5, "float")
                d = unwrap(tokens, 6, "float")
                last_molecule().interactions.append(
                    (system.vsite_3fd, members, [a, d])
                )
            case 3:
                # 3fad
                theta = unwrap(tokens, 5, "float") * math.pi / 180
                d = unwrap(tokens, 6, "float")
                last_molecule().interactions.append(
                    (system.vsite_3fad, members, [theta, d])
                )
            case 4:
                # 3out
                a = unwrap(tokens, 5, "float")
                b = unwrap(tokens, 6, "float")
                c = unwrap(tokens, 7, "float")
                last_molecule().interactions.append(
                    (system.vsite_3out, members, [a, b, c])
                )
            case _:
                raise ValueError(
                    f"Virtual site 3 type {type} not implemented."
                )

    p.add_level("virtual_sites3", process_virtual_sites3)

    def process_virtual_sites4(tokens):
        vid = unwrap(tokens, 0, "index")
        i = unwrap(tokens, 1, "index")
        j = unwrap(tokens, 2, "index")
        k = unwrap(tokens, 3, "index")
        l = unwrap(tokens, 4, "index")
        members = [vid, i, j, k, l]
        type = unwrap(tokens, 5, "int")
        match type:
            case 2:
                a = unwrap(tokens, 6, "float")
                b = unwrap(tokens, 7, "float")
                c = unwrap(tokens, 8, "float")
                last_molecule().interactions.append(
                    (system.vsite_4fdn, members, [a, b, c])
                )

    p.add_level("virtual_sites4", process_virtual_sites4)

    def process_virtual_sitesn(tokens):
        vid = unwrap(tokens, 0, "index")
        members = [vid]
        type = unwrap(tokens, 1, "int")
        match type:
            case 1:
                for c in range(2, len(tokens)):
                    members.append(unwrap(tokens, c, "index"))
                n = len(members) - 1
                weights = [1/n] * n
                last_molecule().interactions.append(
                    (system.vsite_avg, members, weights)
                )
            case 2:
                for c in range(2, len(tokens)):
                    members.append(unwrap(tokens, c, "index"))
                n = len(members) - 1
                last_molecule().interactions.append(
                    (system.vsite_com, members, [])
                )
            case 3:
                weights = []
                wsum = 0.
                if len(tokens) % 2 != 0:
                    raise ValueError("Must have an even number of arguments")
                for c in range(2, len(tokens), 2):
                    index = unwrap(tokens, c, "index")
                    weight = unwrap(tokens, c+1, "float")
                    members.append(index)
                    weights.append(weight)
                    wsum += weight
                weights = [w/wsum for w in weights]
                last_molecule().interactions.append(
                    (system.vsite_avg, members, weights)
                )
            case _:
                raise ValueError(
                    f"Virtua site n type {type} not implemented."
                )

    p.add_level("virtual_sitesn", process_virtual_sitesn)

    # custom additions: rx and frag
    _frag_name = None
    _atom_list = []

    def process_frag(tokens):
        nonlocal _frag_name, _atom_list
        _atom_list = []
        _frag_name = unwrap(tokens, 0, "word")

    p.add_level("frag", process_frag)

    def process_fragatoms(tokens):
        nonlocal _atom_list, _frag_name
        if _frag_name is None:
            raise ValueError("Define a [frag] first.")
        type = unwrap(tokens, 0, {"word", "pattern"})
        name = unwrap(tokens, 1, {"word", "pattern"})
        _atom_list.append((type, name))

    p.add_level("frag_atoms", process_fragatoms)

    def process_fragfrom(tokens):
        nonlocal _frag_name, _atom_list
        if _frag_name is None or len(_atom_list) == 0:
            raise ValueError("Define a [frag] and give [frag_atoms] first.")
        mol = unwrap(tokens, 0, "word")
        parent_ids = []
        for i in range(1, len(tokens)):
            parent_ids.append(unwrap(tokens, i, "index"))
        if len(parent_ids) != len(_atom_list):
            raise ValueError("Wrong number of atoms, expect "
                             f"{len(_atom_list)}, got {len(parent_ids)}")
        frag = topology.new_frag_fragment(_frag_name, mol)
        for i in range(len(parent_ids)):
            frag.atoms.append((parent_ids[i], *(_atom_list[i])))

    p.add_level("frag_from", process_fragfrom)

    # FIXME this is terrible, the current callback architecture is really
    # unsuitable for this
    last_reaction: ReactionTemplate = None

    def require_complete_reaction():
        nonlocal last_reaction
        if last_reaction is None:
            return
        if (last_reaction.r1 is not None and last_reaction.r2 is not None and
                last_reaction.p1 is not None and
                last_reaction.distance_max is not None):
            # complete
            topology.new_reaction(last_reaction)
            last_reaction = None
        else:
            raise ValueError(f"unfinished reaction {last_reaction.name}")

    def process_reaction(tokens):
        nonlocal last_reaction
        if system_defined:
            raise ValueError("[ rx ] must be before [ system ]")
        require_complete_reaction()
        name = unwrap(tokens, 0, "word")
        last_reaction = ReactionTemplate(name, None, None, None, None, [], [])

    p.add_level("rx", process_reaction)

    def process_reactants(tokens):
        nonlocal last_reaction
        r1 = unwrap(tokens, 0, "word")
        r2 = unwrap(tokens, 1, "word")
        last_reaction.r1 = r1
        last_reaction.r2 = r2

    p.add_level("rx_reactants", process_reactants)

    def process_products(tokens):
        nonlocal last_reaction
        p1 = unwrap(tokens, 0, "word")
        last_reaction.p1 = p1

    p.add_level("rx_products", process_products)

    def process_conditions(tokens):
        key, value = unwrap(tokens, 0, "pair")
        match key:
            case "distance_max":
                last_reaction.distance_max = value
            case _:
                raise ValueError(f"Unknown key {key} in [rx_conditions]")

    p.add_level("rx_conditions", process_conditions)

    def process_system(tokens):
        nonlocal system_defined
        require_complete_reaction()
        system_defined = True  # hack so that reactions are complete

    p.add_level("system", process_system)

    # run parser
    ok = p.parse(file, include_dir, defines)
    if not ok:
        raise ValueError("Parsing error")

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
