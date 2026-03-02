"""
A more modular
parser that constructs S* and T* rather than openmm's internal objects
"""

from .parser import Parser, unwrap, TokenParseError, ParseError
from .topstar import TopStar
from .molecule import Molecule
from .reaction_template import ReactionTemplate
from .sysstar import SysStar
from .graph import Graph, GraphAtomType
import math
import logging
import difflib


def DaemonTopFile(
    file, nonbonded, include_dir=None, defines: dict[str, str] = {},
    nlist_cutoff: float = 1.1,
    max_absolute_rate: float | None = None,
    rate_highest_probability: float = 1.0,
    logger=None, respos=None,
    experimental=False
) -> tuple[SysStar, TopStar]:
    """
    Parses a Martini Top file for Gromacs and generates T*, sys and top
    from it. Also parses .frag and .rx files included in the .top file.
    """

    if logger is None:
        logging.basicConfig(filename="out.log", level=logging.INFO)
        logger = logging.getLogger(__name__)
    # field init
    system = SysStar(logger, nonbonded)
    topology = TopStar(
        system, logger, nlist_cutoff,
        max_absolute_rate, rate_highest_probability,
        respos
    )

    # so .top files can be compatible with gromacs
    defines["DAEMON"] = "1"

    # make parser
    p = Parser()

    # state of the parser
    system_defined = False  # whether the [system] happened yet
    last_molecule: Molecule | None = None
    last_graph: Graph | None = None
    last_reaction: ReactionTemplate | None = None

    def this_is_experimental(tokens, index, feature):
        if not experimental:
            raise TokenParseError(
                tokens[index],
                f"{feature} are an experimental feature, it may not work correctly."
                " To enable it, pass the argument experimental=True."
            )

    # a list of parsers for each directive
    def parse_pair(tokens, id, filter) -> int | tuple[int, int]:
        nonlocal last_reaction
        pair = unwrap(tokens, id, filter)
        if type(pair) is int:
            return pair
        if last_reaction is None:
            raise TokenParseError(
                tokens[id],
                "Attempt to pair-index before a [reaction] was defined"
            )
        n_reac = len(last_reaction.reactants)
        idi, namei = pair
        if idi < 0 or idi >= n_reac:
            raise TokenParseError(
                tokens[id],
                f"Reactant index {idi+1} out of range: 1 to {n_reac}."
            )
        graph_name = last_reaction.reactants[idi]
        graph = topology.graphs.get(graph_name)
        if graph is None:
            raise TokenParseError(
                tokens[id],
                f"Reactant graph {graph_name} not found."
            )
        atomi = graph.atom_name_to_index.get(namei)
        if atomi is None:
            raise TokenParseError(
                tokens[id],
                f"Atom {namei} not found in reactant {graph_name}."
            )
        if graph.atoms[atomi][3] not in {
            GraphAtomType.NORMAL, GraphAtomType.OPT
        }:
            raise TokenParseError(
                tokens[id],
                "Attempt to reference a forbidden atom!"
            )
        return idi, atomi

    def process_defaults(tokens):
        nb_type = unwrap(tokens, 0, "int")
        if nb_type != 1:
            raise TokenParseError(
                tokens[0],
                f"Unsupported nonbonded type {nb_type}."
            )
        combination_rule = unwrap(tokens, 1, "int")
        if combination_rule != 2:
            raise TokenParseError(
                tokens[1],
                f"Unsupported combination rule {combination_rule}."
            )
        if len(tokens) > 2:
            raise TokenParseError(
                tokens[2],
                "Too many tokens in [defaults] directive, "
                f"got {len(tokens)}, expect 2."
            )

    p.add_level("defaults", process_defaults)

    def process_moltype(tokens):
        nonlocal last_molecule
        name = unwrap(tokens, 0, "word")
        nrexcl = unwrap(tokens, 1, "int")
        if nrexcl != 1:
            raise TokenParseError(
                tokens[1],
                f"nrexcl is {nrexcl}, only nrexcl=1 is supported."
            )
        last_molecule = topology.new_molecule(name)

    p.add_level("moleculetype", process_moltype)

    def molecules_start():
        if not system_defined:
            raise ParseError(
                "[molecules] must come after [system]."
            )

    def process_molecule(tokens):
        name = unwrap(tokens, 0, "word")
        if topology.molecules.get(name) is None:
            raise TokenParseError(
                tokens[0],
                f"Undefined molecule type {name}."
            )
        count = unwrap(tokens, 1, "int")
        for i in range(count):
            topology.instantiate(name)

    p.add_level("molecules", process_molecule, start=molecules_start)

    def assert_last_molecule():
        nonlocal last_molecule
        if last_molecule is None:
            raise ParseError("No [moleculetype] given.")

    def process_atoms(tokens):
        nonlocal last_molecule
        if last_molecule.index_type != "index":
            raise ParseError("[atoms] only valid in [moleculetype]")
        id = unwrap(tokens, 0, "index")
        type = unwrap(tokens, 1, "word")
        resnum = unwrap(tokens, 2, "int")
        resname = unwrap(tokens, 3, "word")
        atomname = unwrap(tokens, 4, "word")
        charge_group_num = unwrap(tokens, 5, "int")
        charge = unwrap(tokens, 6, "float") if len(tokens) > 6 else None
        mass = unwrap(tokens, 7, "float") if len(tokens) > 7 else None
        atom_index = len(last_molecule.atoms)
        if id != atom_index:
            raise TokenParseError(
                tokens[0],
                "Bad atom ID, are they out of order?"
                f" got id {id} but expected {atom_index}"
            )
        last_molecule.atoms.append(
            (type, resnum, resname, atomname, charge_group_num, charge, mass)
        )

    p.add_level("atoms", process_atoms, start=assert_last_molecule)

    def process_bonds(tokens):
        nonlocal last_molecule
        index_type = last_molecule.index_type
        i = parse_pair(tokens, 0, index_type)
        j = parse_pair(tokens, 1, index_type)
        type = unwrap(tokens, 2, "int")
        length = None
        if type != 5:
            # type 5 takes no length
            length = unwrap(tokens, 3, "float")
        if type == 1 or type == 6:
            # harmonic bond / harmonic potential
            force = unwrap(tokens, 4, "float")
            last_molecule.interactions.append(
                (system.harmonic_bond, [i, j], [length, force])
            )
        elif type == 2:
            # G96 bond
            force = unwrap(tokens, 4, "float")
            last_molecule.interactions.append(
                (system.g96_bond, [i, j], [length, force])
            )
        elif type == 3:
            # morse
            D = unwrap(tokens, 4, "float")
            beta = unwrap(tokens, 5, "float")
            last_molecule.interactions.append(
                (system.morse_bond, [i, j], [length, D, beta])
            )
        elif type == 4:
            # cubic bond
            kb = unwrap(tokens, 4, "float")
            kcub = unwrap(tokens, 5, "float")
            last_molecule.interactions.append(
                (system.cubic_bond, [i, j], [length, kb, kcub])
            )
        elif type == 5:
            # connection - only generate exclusions and add a dummy interaction
            last_molecule.interactions.append(
                (system.connection, [i, j], [])
            )
        elif type == 7:
            # FENE (finitely extensible nonlinear elastic) bond
            force = unwrap(tokens, 4, "float")
            last_molecule.interactions.append(
                (system.fene_bond, [i, j], [length, force])
            )
        elif type == 10:
            # restraint potential
            low = length
            up1 = unwrap(tokens, 4, "float")
            up2 = unwrap(tokens, 5, "float")
            k = unwrap(tokens, 6, "float")
            last_molecule.interactions.append(
                (system.distance_restraint, [i, j], [low, up1, up2, k])
            )
        else:
            raise TokenParseError(
                tokens[2],
                f"Unsupported bond function type {type}."
            )
        if type != 6:
            # type 6 does not generate exclusions
            # nrexcl other than 1 is not supported
            last_molecule.add_exclusion(i, j)

    p.add_level("bonds", process_bonds, start=assert_last_molecule)

    def process_angles(tokens):
        nonlocal last_molecule
        index_type = last_molecule.index_type
        i = parse_pair(tokens, 0, index_type)
        j = parse_pair(tokens, 1, index_type)
        k = parse_pair(tokens, 2, index_type)
        type = unwrap(tokens, 3, "int")
        if type == 1:
            # Harmonic angle
            theta = unwrap(tokens, 4, "degree")
            force = unwrap(tokens, 5, "float")
            last_molecule.interactions.append(
                (system.harmonic_angle, [i, j, k], [theta, force])
            )
        elif type == 2:
            # G96 angle
            theta = unwrap(tokens, 4, "degree")
            force = unwrap(tokens, 5, "float")
            last_molecule.interactions.append(
                (system.g96_angle, [i, j, k], [theta, force])
            )
        elif type == 3:
            # cross bond bond
            r1 = unwrap(tokens, 4, "float")
            r2 = unwrap(tokens, 5, "float")
            force = unwrap(tokens, 6, "float")
            last_molecule.interactions.append(
                (system.cross_bond_bond, [i, j, k], [r1, r2, force])
            )
        elif type == 4:
            # cross bond angle
            r1 = unwrap(tokens, 4, "float")
            r2 = unwrap(tokens, 5, "float")
            r3 = unwrap(tokens, 6, "float")
            force = unwrap(tokens, 7, "float")
            last_molecule.interactions.append(
                (system.cross_bond_angle, [i, j, k], [r1, r2, r3, force])
            )
        elif type == 5:
            # Urey-Bradley
            theta = unwrap(tokens, 4, "degree")
            force = unwrap(tokens, 5, "float")
            r13 = unwrap(tokens, 6, "float")
            k_UB = unwrap(tokens, 7, "float")
            last_molecule.interactions.append(
                (system.urey_bradley, [i, j, k], [theta, force, r13, k_UB])
            )
        elif type == 6:
            # Quartic angle
            # theta
            params = [unwrap(tokens, 4, "degree")]
            # C0 to C4
            for x in range(5):
                params.append(unwrap(tokens, 5+x, "float"))
            last_molecule.interactions.append(
                (system.quartic_angle, [i, j, k], params)
            )
        elif type == 9:
            # Linear angle
            a = unwrap(tokens, 4, "float")
            force = unwrap(tokens, 5, "float")
            last_molecule.interactions.append(
                (system.linear_angle, [i, j, k], [a, force])
            )
        elif type == 10:
            # Restricted angle
            theta = unwrap(tokens, 4, "degree")
            force = unwrap(tokens, 5, "float")
            last_molecule.interactions.append(
                (system.restricted_angle, [i, j, k], [theta, force])
            )
        else:
            raise TokenParseError(
                tokens[3],
                f"Unsupported angle type {type}."
            )

    p.add_level("angles", process_angles, start=assert_last_molecule)

    def process_dihedrals(tokens):
        nonlocal last_molecule
        index_type = last_molecule.index_type
        i = parse_pair(tokens, 0, index_type)
        j = parse_pair(tokens, 1, index_type)
        k = parse_pair(tokens, 2, index_type)
        L = parse_pair(tokens, 3, index_type)
        type = unwrap(tokens, 4, "int")
        match type:
            case 1 | 4 | 9:
                # proper dihedral | proper dihedral (multiple)
                theta = unwrap(tokens, 5, "degree")
                force = unwrap(tokens, 6, "float")
                multiplicity = unwrap(tokens, 7, "int", 1)
                last_molecule.interactions.append(
                    (system.proper_dihedral, [i, j, k, L],
                     [theta, force, multiplicity])
                )
            case 2:
                # improper
                theta = unwrap(tokens, 5, "degree")
                force = unwrap(tokens, 6, "float")
                theta = theta - math.tau if theta > math.pi else theta
                last_molecule.interactions.append(
                    (system.improper_dihedral, [i, j, k, L],
                     [theta, force])
                )
            case 3:
                # ryckaert-bellemans = RB
                params = []
                for c in range(6):
                    # C0 to C5
                    params.append(unwrap(tokens, 5+c, "float"))
                last_molecule.interactions.append(
                    (system.rb_torsion, [i, j, k, L], params)
                )
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
                last_molecule.interactions.append(
                    (system.rb_torsion, [i, j, k, L], rb_params)
                )
            case 10:
                # restricted dihedral
                theta = unwrap(tokens, 5, "degree")
                force = unwrap(tokens, 6, "float")
                last_molecule.interactions.append(
                    (system.restricted_dihedral, [i, j, k, L], [theta, force])
                )
            case 11:
                # combined bending-torsion potential
                force = unwrap(tokens, 5, "float")
                params = [force]
                for c in range(5):
                    # a0 to a4
                    params.append(unwrap(tokens, 6+c, "float"))
                last_molecule.interactions.append(
                    (system.combined_bending_torsion, [i, j, k, L], params)
                )
            case 101:
                this_is_experimental(
                    tokens, 4, "Periodic gaussian dihedrals"
                )
                # custom type - periodic gaussian
                theta = unwrap(tokens, 5, "degree")
                depth = unwrap(tokens, 6, "float")
                force = unwrap(tokens, 7, "float")
                last_molecule.interactions.append((
                    system.periodic_gaussian, [i, j, k, L],
                    [theta, depth, force]
                ))
            case _:
                raise TokenParseError(
                    tokens[4],
                    f"Unsupported dihedral function type {type}."
                )

    p.add_level("dihedrals", process_dihedrals, start=assert_last_molecule)

    def process_exclusions(tokens):
        nonlocal last_molecule
        index_type = last_molecule.index_type
        i = parse_pair(tokens, 0, index_type)
        j = parse_pair(tokens, 1, index_type)
        last_molecule.add_exclusion(i, j)
        for k in range(2, len(tokens)):
            c = unwrap(tokens, k, index_type)
            last_molecule.add_exclusion(i, c)

    p.add_level("exclusions", process_exclusions, start=assert_last_molecule)

    def process_constraints(tokens):
        nonlocal last_molecule
        index_type = last_molecule.index_type
        i = parse_pair(tokens, 0, index_type)
        j = parse_pair(tokens, 1, index_type)
        type = unwrap(tokens, 2, "int")
        length = unwrap(tokens, 3, "float")
        if type == 1 or type == 2:
            last_molecule.interactions.append(
                (system.constraint, [i, j], [length])
            )
            if type == 1:
                # type 2 doesn't generate exclusions
                last_molecule.add_exclusion(i, j)
        else:
            raise TokenParseError(
                tokens[2],
                f"Unsupported constraint type {type}."
            )

    p.add_level("constraints", process_constraints, start=assert_last_molecule)

    def unwrap_atomtype(tokens, index):
        type = unwrap(tokens, index, "word")
        if system._atom_types.get(type) is None:
            raise TokenParseError(
                tokens[index],
                f"Unknown atom type {type}."
            )
        return type

    def process_pairtypes(tokens):
        this_is_experimental(
            tokens, 0, "Pairs"
        )
        t1 = unwrap_atomtype(tokens, 0)
        t2 = unwrap_atomtype(tokens, 1)
        type = unwrap(tokens, 2, "int")
        if type != 1:
            raise TokenParseError(
                tokens[2],
                f"Unsupported pairs type {type}."
            )
        sigma = unwrap(tokens, 3, "float")
        epsilon = unwrap(tokens, 4, "float")
        system.pairs.add_type(t1, t2, sigma, epsilon)

    p.add_level("pairtypes", process_pairtypes)

    def process_pairs(tokens):
        nonlocal last_molecule
        index_type = last_molecule.index_type
        i = parse_pair(tokens, 0, index_type)
        j = parse_pair(tokens, 1, index_type)
        type = unwrap(tokens, 2, "int")
        if type != 1:
            raise TokenParseError(
                tokens[2],
                f"Unsupported pairs type {type}."
            )
        if len(tokens) >= 5:
            sigma = unwrap(tokens, 3, "float")
            epsilon = unwrap(tokens, 4, "float")
            last_molecule.interactions.append(
                (system.pairs, [i, j], [sigma, epsilon])
            )
        else:
            last_molecule.interactions.append(
                (system.pairs, [i, j], [])
            )

    p.add_level("pairs", process_pairs, start=assert_last_molecule)

    def process_cmaptypes(tokens):
        parts = [
            unwrap_atomtype(tokens, i) for i in range(5)
        ]
        type = unwrap(tokens, 5, "int")
        if type != 1:
            raise TokenParseError(
                tokens[5],
                f"Unsupported cmap type {type}."
            )
        size = [
            unwrap(tokens, 6, "int"),
            unwrap(tokens, 7, "int")
        ]
        if size[0] != size[1]:
            raise TokenParseError(
                tokens[7],
                "Non-square CMAPs are not supported."
            )
        if size < 8:
            logger.warn(
                f"CMAPs of size {size}x{size} are small. This might result in "
                "slightly different behavior between GROMACS and OpenMM."
            )
        params = [
            unwrap(tokens, 8+i, "float")
            for i in range(size[0] * size[1])
        ]
        system.cmap.add_type(parts, [*size, *params])

    p.add_level("cmaptypes", process_cmaptypes)

    def process_cmap(tokens):
        this_is_experimental(
            tokens, 0, "Cmaps"
        )
        nonlocal last_molecule
        index_type = last_molecule.index_type
        members = [
            parse_pair(tokens, i, index_type)
            for i in range(5)
        ]
        type = unwrap(tokens, 5, "int")
        if type != 1:
            raise TokenParseError(
                tokens[5],
                f"Unsupported cmap type {type}."
            )
        # GromacsTopFile in openmm does support cmap types
        # here, but GROMACS seems to require [cmaptypes]
        # so we only support [cmaptypes]
        last_molecule.interactions.append((
            system.cmap, members, []
        ))

    p.add_level("cmap", process_cmap, start=assert_last_molecule)

    def process_reactive_type(tokens):
        this_is_experimental(
            tokens, 0, "Custom reactive types"
        )
        filter = unwrap(tokens, 0, "word")

        if system.custom_reactive.reactive_types.get(filter) is not None:
            raise TokenParseError(
                tokens[0],
                f"Reactive type {filter} already defined."
            )

        i = 1

        """
        sigma = 0.
        epsilon = 0.
        barrier_height = 0.
        barrier_position = 0.
        barrier_width = 0.01
        well_position = 0.
        well_depth = 0.
        well_width = 0.01
        repulsive_strength = 0.
        angle = 0.
        angle_strength = 0.
        dihedral = 0.
        dihedral_depth = 0.
        dihedral_width = 0.01
        """

        keys = set()
        while i < len(tokens):
            key = unwrap(tokens, i, "word")
            if key in keys:
                raise TokenParseError(
                    tokens[i],
                    f"Duplicate entry for key {key}."
                )
            keys.add(key)
            match key:
                case "LJ_negate":
                    sigma = unwrap(tokens, i+1, "float")
                    epsilon = unwrap(tokens, i+2, "float")
                    i += 3
                case "barrier":
                    barrier_position = unwrap(tokens, i+2, "float")
                    barrier_height = unwrap(tokens, i+1, "float")
                    barrier_width = unwrap(tokens, i+3, "float")
                    i += 4
                case "well":
                    well_position = unwrap(tokens, i+1, "float")
                    well_depth = unwrap(tokens, i+2, "float")
                    well_width = unwrap(tokens, i+3, "float")
                    repulsive_strength = unwrap(tokens, i+4, "float")
                    i += 5
                case "angle":
                    angle = unwrap(tokens, i+1, "degree")
                    angle_strength = unwrap(tokens, i+2, "float")
                    i += 3
                case "dihedral":
                    dihedral = unwrap(tokens, i+1, "degree")
                    if dihedral < 0. or dihedral > math.pi:
                        raise TokenParseError(
                            tokens[i+1],
                            "Value must be between 0 and 180 degrees."
                        )
                    dihedral_depth = unwrap(tokens, i+2, "float")
                    dihedral_width = unwrap(tokens, i+3, "float")
                    i += 4
                case _:
                    raise TokenParseError(
                        tokens[i],
                        f"Unknown key {key}."
                    )

        # TODO defaults, so it's not all mandatory
        mandatory = {"LJ_negate", "barrier", "well", "angle", "dihedral"}
        remaining = mandatory - keys
        if len(remaining) > 0:
            raise ParseError(
                f"Missing entries {remaining}."
            )
        system.custom_reactive.reactive_types[filter] = (
            sigma, epsilon,
            barrier_position, barrier_height, barrier_width,
            well_position, well_depth, well_width, repulsive_strength,
            angle, angle_strength,
            dihedral, dihedral_depth, dihedral_width
        )

    p.add_level(
        "reactive_types", process_reactive_type
    )

    def process_custom_reactive(tokens):
        this_is_experimental(
            tokens, 0, "Custom reactive types"
        )
        nonlocal last_molecule
        # just a test for now, hardcoded constants
        index_type = last_molecule.index_type
        i = parse_pair(tokens, 0, index_type)
        j = parse_pair(tokens, 1, index_type)
        k = parse_pair(tokens, 2, index_type)

        type = unwrap(tokens, 3, "word")
        if type not in {"d", "a", "s"}:
            raise TokenParseError(
                tokens[3],
                f"Must be one of {'d', 'a', 's'}, got {type}."
                " d=donor, a=acceptor, s=symmetric."
            )
        filter = unwrap(tokens, 4, "word")
        if system.custom_reactive.reactive_types.get(filter) is None:
            raise TokenParseError(
                tokens[4],
                f"Undefined reactive type {filter}."
                "Define it in [reactive_types] first."
            )

        last_molecule.interactions.append((
            system.custom_reactive, [i, j, k], [type, filter]
        ))

    p.add_level(
        "reactive_group", process_custom_reactive,
        start=assert_last_molecule
    )

    def process_atomtypes(tokens):
        if len(tokens) != 6:
            raise ParseError(
                "Only atomtypes lines formatted as type, m, q,"
                "atom_type, sigma, epsilon are supported."
            )
        type = unwrap(tokens, 0, "word")
        mass = unwrap(tokens, 1, "float")
        charge = unwrap(tokens, 2, "float")
        atom_type = unwrap(tokens, 3, "word")
        if atom_type != "A":
            raise TokenParseError(
                tokens[3],
                "Only 'A' atom type supported."
            )
        sigma = unwrap(tokens, 4, "float")
        epsilon = unwrap(tokens, 5, "float")
        if sigma != 0.0:
            raise TokenParseError(
                tokens[4],
                "Only sigma=0 is supported in [atomtypes]."
            )
        if epsilon != 0.0:
            raise TokenParseError(
                tokens[5],
                "Only epsilon=0 is supported in [atomtypes]."
            )
        system.add_atom_type(type, charge, mass)

    p.add_level("atomtypes", process_atomtypes)

    def process_position_restraints(tokens):
        nonlocal last_molecule
        if len(tokens) != 5:
            raise ParseError(
                "Only position_restraint lines formatted as "
                "<index> <type> <kx> <ky> <kz> are supported."
            )
        type = unwrap(tokens, 1, "int")
        if type != 1:
            raise TokenParseError(
                tokens[1],
                f"Unsupported position restraint type {type}."
            )
        id = unwrap(tokens, 0, "index")
        kx = unwrap(tokens, 2, "float")
        ky = unwrap(tokens, 3, "float")
        kz = unwrap(tokens, 4, "float")
        # x0, y0, z0 can only be set during instantiation - it is part of T*
        last_molecule.posres.append((id, kx, ky, kz))

    p.add_level(
        "position_restraints", process_position_restraints,
        start=assert_last_molecule
    )

    def process_nonbond_params(tokens):
        type1 = unwrap(tokens, 0, "word")
        type2 = unwrap(tokens, 1, "word")
        funct = unwrap(tokens, 2, "int")
        if funct != 1:
            raise TokenParseError(
                tokens[2],
                f"Unsupported function type {type}."
            )
        sigma = unwrap(tokens, 3, "float")
        epsilon = unwrap(tokens, 4, "float")
        system.add_nb_type(type1, type2, sigma, epsilon)

    p.add_level("nonbond_params", process_nonbond_params)

    def process_virtual_sites1(tokens):
        nonlocal last_molecule
        vid = parse_pair(tokens, 0, last_molecule.index_type)
        member = parse_pair(tokens, 1, last_molecule.index_type)
        type = unwrap(tokens, 2, "int")
        if type != 1:
            raise TokenParseError(
                tokens[2],
                f"Virtual site 1 type {type} not supported."
            )
        last_molecule.interactions.append(
            (system.vsite1, [vid, member], [])
        )
        # OpenMM does not like overlapping atoms otherwise
        last_molecule.add_exclusion(vid, member)

    p.add_level(
        "virtual_sites1", process_virtual_sites1, start=assert_last_molecule
    )

    def process_virtual_sites2(tokens):
        nonlocal last_molecule
        index_type = last_molecule.index_type
        vid = parse_pair(tokens, 0, index_type)
        i = parse_pair(tokens, 1, index_type)
        j = parse_pair(tokens, 2, index_type)
        members = [vid, i, j]
        type = unwrap(tokens, 3, "int")
        match type:
            case 1:
                a = unwrap(tokens, 4, "float")
                weights = [1-a, a]
                last_molecule.interactions.append(
                    (system.vsite2, members, weights)
                )
            case 2:
                d = unwrap(tokens, 4, "float")
                last_molecule.interactions.append(
                    (system.vsite_2fd, members, [d])
                )
            case _:
                raise TokenParseError(
                    tokens[3],
                    f"Virtual site 2 type {type} unsupported."
                )

    p.add_level(
        "virtual_sites2", process_virtual_sites2,
        start=assert_last_molecule
    )

    def process_virtual_sites3(tokens):
        nonlocal last_molecule
        index_type = last_molecule.index_type
        vid = parse_pair(tokens, 0, index_type)
        i = parse_pair(tokens, 1, index_type)
        j = parse_pair(tokens, 2, index_type)
        k = parse_pair(tokens, 3, index_type)
        members = [vid, i, j, k]
        type = unwrap(tokens, 4, "int")
        match type:
            case 1:
                # 3
                a = unwrap(tokens, 5, "float")
                b = unwrap(tokens, 6, "float")
                weights = [1 - a - b, a, b]
                last_molecule.interactions.append(
                    (system.vsite3, members, weights)
                )
            case 2:
                # 3fd
                a = unwrap(tokens, 5, "float")
                d = unwrap(tokens, 6, "float")
                last_molecule.interactions.append(
                    (system.vsite_3fd, members, [a, d])
                )
            case 3:
                # 3fad
                theta = unwrap(tokens, 5, "degree")
                d = unwrap(tokens, 6, "float")
                last_molecule.interactions.append(
                    (system.vsite_3fad, members, [theta, d])
                )
            case 4:
                # 3out
                a = unwrap(tokens, 5, "float")
                b = unwrap(tokens, 6, "float")
                c = unwrap(tokens, 7, "float")
                last_molecule.interactions.append(
                    (system.vsite_3out, members, [a, b, c])
                )
            case _:
                raise TokenParseError(
                    tokens[4],
                    f"Virtual site 3 type {type} unsupported."
                )

    p.add_level(
        "virtual_sites3", process_virtual_sites3,
        start=assert_last_molecule
    )

    def process_virtual_sites4(tokens):
        nonlocal last_molecule
        index_type = last_molecule.index_type
        vid = parse_pair(tokens, 0, index_type)
        i = parse_pair(tokens, 1, index_type)
        j = parse_pair(tokens, 2, index_type)
        k = parse_pair(tokens, 3, index_type)
        L = parse_pair(tokens, 4, index_type)
        members = [vid, i, j, k, L]
        type = unwrap(tokens, 5, "int")
        match type:
            case 2:
                a = unwrap(tokens, 6, "float")
                b = unwrap(tokens, 7, "float")
                c = unwrap(tokens, 8, "float")
                last_molecule.interactions.append(
                    (system.vsite_4fdn, members, [a, b, c])
                )
            case _:
                raise TokenParseError(
                    tokens[5],
                    f"Virtual site 4 type {type} unsupported."
                )

    p.add_level(
        "virtual_sites4", process_virtual_sites4,
        start=assert_last_molecule
    )

    def process_virtual_sitesn(tokens):
        nonlocal last_molecule
        index_type = last_molecule.index_type
        vid = parse_pair(tokens, 0, index_type)
        members = [vid]
        type = unwrap(tokens, 1, "int")
        match type:
            case 1:
                for c in range(2, len(tokens)):
                    members.append(parse_pair(tokens, c, index_type))
                n = len(members) - 1
                weights = [1/n] * n
                last_molecule.interactions.append(
                    (system.vsite_avg, members, weights)
                )
            case 2:
                for c in range(2, len(tokens)):
                    members.append(parse_pair(tokens, c, index_type))
                n = len(members) - 1
                last_molecule.interactions.append(
                    (system.vsite_com, members, [])
                )
            case 3:
                weights = []
                wsum = 0.
                if len(tokens) % 2 != 0:
                    raise ParseError("Must have an even number of arguments")
                for c in range(2, len(tokens), 2):
                    index = parse_pair(tokens, c, index_type)
                    weight = unwrap(tokens, c+1, "float")
                    members.append(index)
                    weights.append(weight)
                    wsum += weight
                weights = [w/wsum for w in weights]
                last_molecule.interactions.append(
                    (system.vsite_avg, members, weights)
                )
            case _:
                raise TokenParseError(
                    tokens[1],
                    f"Virtual site n type {type} unsupported."
                )
        # OpenMM does not like overlapping particles without exclusions
        if len(members) == 2:
            last_molecule.add_exclusion(members[0], members[1])

    p.add_level(
        "virtual_sitesn", process_virtual_sitesn,
        start=assert_last_molecule
    )

    def graph_start():
        nonlocal last_graph
        last_graph = Graph(None)

    def process_graph(tokens):
        nonlocal last_graph
        keyword = unwrap(tokens, 0, "word")
        if keyword == "name":
            if last_graph.name is not None:
                raise TokenParseError(
                    tokens[0],
                    f"The graph already has a name {last_graph.name}."
                )
            last_graph.name = unwrap(tokens, 1, "word")
        elif keyword == "atom":
            part_id = unwrap(tokens, 1, "word")
            name_pat = unwrap(tokens, 2, "pattern")
            type_pat = unwrap(tokens, 3, "pattern")
            type = GraphAtomType.NORMAL
            last_graph.atoms.append((part_id, name_pat, type_pat, type))
        elif keyword == "atom?":
            part_id = unwrap(tokens, 1, "word")
            name_pat = unwrap(tokens, 2, "pattern")
            type_pat = unwrap(tokens, 3, "pattern")
            type = GraphAtomType.OPT
            last_graph.atoms.append((part_id, name_pat, type_pat, type))
        elif keyword == "atom!":
            part_id = unwrap(tokens, 1, "word")
            name_pat = unwrap(tokens, 2, "pattern")
            type_pat = unwrap(tokens, 3, "pattern")
            type = GraphAtomType.NOT
            last_graph.atoms.append((part_id, name_pat, type_pat, type))
        elif keyword == "molecule":
            mols = [unwrap(tokens, n, "word") for n in range(1, len(tokens))]
            last_graph.molecules += mols
        elif keyword == "equivalent":
            parts = set(
                unwrap(tokens, n, "word") for n in range(1, len(tokens))
            )
            last_graph.equivalents.append(parts)
        else:
            if keyword in system.filters:
                parts = [unwrap(tokens, i, "word") for i in range(1, len(tokens))]
                last_graph.interactions.append((keyword, parts))
            else:
                possibilities = system.filters | {
                    "atom", "atom?", "atom!", "molecule", "name"
                }
                close_matches = difflib.get_close_matches(
                    keyword, possibilities, 1
                )
                raise TokenParseError(
                    tokens[0],
                    f"Keyword {keyword} in graph unrecognized." +
                    (
                        f" Perhaps you meant {close_matches[0]}?"
                        if len(close_matches) > 0 else ""
                    )
                )

    def graph_end():
        nonlocal last_graph, topology
        if last_graph.name is None:
            ParseError("Graph name is mandatory, but was not given.")
        last_graph.finish_init()
        topology.new_graph(last_graph)
        last_graph = None

    p.add_level("graph", process_graph, start=graph_start, end=graph_end)
    # alias:
    p.add_level("frag", process_graph, start=graph_start, end=graph_end)

    def require_complete_reaction():
        nonlocal last_reaction
        if last_reaction is None:
            return
        if last_reaction.is_complete():
            # complete
            topology.new_reaction(last_reaction)
            last_reaction = None
        else:
            raise ParseError(f"Unfinished reaction {last_reaction.name}.")

    def process_reaction(tokens):
        nonlocal last_reaction, last_molecule
        if system_defined:
            raise ParseError("[ reaction ] must be before [ system ]")
        require_complete_reaction()
        name = unwrap(tokens, 0, "word")
        last_reaction = ReactionTemplate(name)
        last_molecule = topology.new_molecule(name)
        last_molecule.index_type = "pair"

    p.add_level("reaction", process_reaction)
    # alias:
    p.add_level("rx", process_reaction)

    def process_reactants(tokens):
        nonlocal last_reaction
        rxs = [unwrap(tokens, i, "word") for i in range(len(tokens))]
        for i, rx in enumerate(rxs):
            if topology.graphs.get(rx) is None:
                raise TokenParseError(
                    tokens[i],
                    f"Undefined reactant graph name {rx}."
                    " Please define the graph first."
                )
        if len(tokens) > 4:
            raise ParseError("Only up to 4 reactants are allowed.")
        if len(last_reaction.reactants) > 0:
            raise ParseError("Only one set of reactants per reaction.")
        last_reaction.reactants = rxs

    p.add_level("reactants", process_reactants)
    # alias:
    p.add_level("reactant", process_reactants)

    def process_conditions(tokens):
        nonlocal last_reaction
        n_reac = len(last_reaction.reactants)
        if n_reac == 0:
            raise ParseError("[reactants] must come before [conditions]")
        key = unwrap(tokens, 0, "word")
        angles = {
            (a, b, c, d, e, f)
            for a, b, c, d, e, f, _, _ in last_reaction.angle_limits
        }
        dihedrals = {
            (a, b, c, d, e, f, g, h)
            for a, b, c, d, e, f, g, h, _, _ in last_reaction.dihedral_limits
        }
        match key:
            case "r_max":
                idi, atomi = parse_pair(tokens, 1, "pair")
                idj, atomj = parse_pair(tokens, 2, "pair")
                cutoff = unwrap(tokens, 3, "positive")
                last_reaction.distance_max.append(
                    (idi, atomi, idj, atomj, cutoff)
                )
            case "r_min":
                idi, atomi = parse_pair(tokens, 1, "pair")
                idj, atomj = parse_pair(tokens, 2, "pair")
                cutoff = unwrap(tokens, 3, "positive")
                last_reaction.distance_min.append(
                    (idi, atomi, idj, atomj, cutoff)
                )
            case "angle":
                idi, atomi = parse_pair(tokens, 1, "pair")
                idj, atomj = parse_pair(tokens, 2, "pair")
                idk, atomk = parse_pair(tokens, 3, "pair")
                if (idi, atomi, idj, atomj, idk, atomk) in angles:
                    raise ParseError(
                        "There is already angle conditions defined for"
                        " the atoms:"
                        f" {unwrap(tokens, 1, 'pair')},"
                        f" {unwrap(tokens, 2, 'pair')},"
                        f" {unwrap(tokens, 3, 'pair')}."
                    )
                # TODO finish this
            case "angle_not":
                idi, atomi = parse_pair(tokens, 1, "pair")
                idj, atomj = parse_pair(tokens, 2, "pair")
                idk, atomk = parse_pair(tokens, 3, "pair")
                theta_min = unwrap(tokens, 4, "degree")
                theta_max = unwrap(tokens, 5, "degree")
                last_reaction.angle_limits.append((
                    idi, atomi, idj, atomj, idk, atomk,
                    math.cos(theta_min), math.cos(theta_max)
                ))
            case "angle_min":
                idi, atomi = parse_pair(tokens, 1, "pair")
                idj, atomj = parse_pair(tokens, 2, "pair")
                idk, atomk = parse_pair(tokens, 3, "pair")
                theta_min = unwrap(tokens, 4, "degree")
                last_reaction.angle_limits.append((
                    idi, atomi, idj, atomj, idk, atomk,
                    math.cos(theta_min), math.cos(math.pi)
                ))
            case "angle_max":
                idi, atomi = parse_pair(tokens, 1, "pair")
                idj, atomj = parse_pair(tokens, 2, "pair")
                idk, atomk = parse_pair(tokens, 3, "pair")
                theta_max = unwrap(tokens, 4, "degree")
                last_reaction.angle_limits.append((
                    idi, atomi, idj, atomj, idk, atomk,
                    math.cos(0.), math.cos(theta_max)
                ))
            case "angle_between":
                idi, atomi = parse_pair(tokens, 1, "pair")
                idj, atomj = parse_pair(tokens, 2, "pair")
                idk, atomk = parse_pair(tokens, 3, "pair")
                theta_min = unwrap(tokens, 4, "degree")
                theta_max = unwrap(tokens, 5, "degree")
                last_reaction.angle_limits.append((
                    idi, atomi, idj, atomj, idk, atomk,
                    math.cos(0.), math.cos(theta_min)
                ))
                last_reaction.angle_limits.append((
                    idi, atomi, idj, atomj, idk, atomk,
                    math.cos(theta_max), math.cos(math.pi)
                ))
            case "dihedral_not":
                idi, atomi = parse_pair(tokens, 1, "pair")
                idj, atomj = parse_pair(tokens, 2, "pair")
                idk, atomk = parse_pair(tokens, 3, "pair")
                idl, atoml = parse_pair(tokens, 4, "pair")
                theta_min = unwrap(tokens, 5, "degree")
                theta_max = unwrap(tokens, 6, "degree")
                # periodicity consideration
                if theta_max > theta_min:
                    last_reaction.dihedral_limits.append((
                        idi, atomi, idj, atomj, idk, atomk, idl, atoml,
                        theta_min, theta_max
                    ))
                else:
                    # minus one / add one for inclusive ranges
                    last_reaction.dihedral_limits.append((
                        idi, atomi, idj, atomj, idk, atomk, idl, atoml,
                        0., theta_max
                    ))
                    last_reaction.dihedral_limits.append((
                        idi, atomi, idj, atomj, idk, atomk, idl, atoml,
                        theta_min, math.tau
                    ))
            case "dihedral_between":
                idi, atomi = parse_pair(tokens, 1, "pair")
                idj, atomj = parse_pair(tokens, 2, "pair")
                idk, atomk = parse_pair(tokens, 3, "pair")
                idl, atoml = parse_pair(tokens, 4, "pair")
                theta_max = unwrap(tokens, 5, "degree")
                theta_min = unwrap(tokens, 6, "degree")
                # periodicity consideration
                if theta_max > theta_min:
                    last_reaction.dihedral_limits.append((
                        idi, atomi, idj, atomj, idk, atomk, idl, atoml,
                        theta_min, theta_max
                    ))
                else:
                    # minus one / add one for inclusive ranges
                    last_reaction.dihedral_limits.append((
                        idi, atomi, idj, atomj, idk, atomk, idl, atoml,
                        0., theta_max
                    ))
                    last_reaction.dihedral_limits.append((
                        idi, atomi, idj, atomj, idk, atomk, idl, atoml,
                        theta_min, math.tau
                    ))
            case "rate":
                this_is_experimental(
                    tokens, 0, "Relative rate control"
                )
                last_reaction.relative_rate = unwrap(tokens, 1, "positive")
            case _:
                raise TokenParseError(
                    tokens[0],
                    f"Unknown reaction condition {key}."
                )

    p.add_level("conditions", process_conditions)
    # alias:
    p.add_level("condition", process_conditions)

    def process_break(tokens):
        nonlocal last_reaction
        n_reac = len(last_reaction.reactants)
        if n_reac == 0:
            raise ParseError("[reactants] must come before [break]")

        atoms = []
        for i in range(len(tokens)):
            atoms.append(parse_pair(tokens, i, "pair"))
        last_reaction.break_groups.append(atoms)

    p.add_level("break", process_break)

    def process_update(tokens):
        nonlocal last_reaction
        n_reac = len(last_reaction.reactants)
        if n_reac == 0:
            raise ParseError("[reactants] must come before [update]")

        atoms = []
        for i in range(len(tokens)):
            atoms.append(parse_pair(tokens, i, "pair"))
        last_reaction.update_groups.append(atoms)

    p.add_level("update", process_update)

    def process_redefine(tokens):
        nonlocal last_reaction
        n_reac = len(last_reaction.reactants)
        if n_reac == 0:
            raise ParseError("[reactants] must come before [redefine]")
        id, atomid = parse_pair(tokens, 0, "pair")
        changes = set()
        if len(tokens) % 2 != 1:
            raise ParseError(
                "[redefine] must contain an atom, and a list of pairs of"
                " properties and new values. Found an even number of"
                " tokens, expected an odd number."
            )
        for i in range((len(tokens)-1) // 2):
            word = unwrap(tokens, 1+i*2, "word")
            if word in changes:
                raise TokenParseError(
                    tokens[1+i*2],
                    f"Duplicate entry {word}."
                )
            changes.add(word)
            match word:
                # TODO unify this into a redefine
                case "name":
                    last_reaction.renames.append((
                        id, atomid, unwrap(tokens, 2+i*2, "word")
                    ))
                case "type":
                    last_reaction.retypes.append((
                        id, atomid, unwrap(tokens, 2+i*2, "word")
                    ))
                case "charge":
                    last_reaction.recharges.append((
                        id, atomid, unwrap(tokens, 2+i*2, "float")
                    ))
                case "mass":
                    last_reaction.remasses.append((
                        id, atomid, unwrap(tokens, 2+i*2, "float")
                    ))
                case _:
                    raise TokenParseError(
                        tokens[1+i*2],
                        f"Unknown atom property {word}."
                    )

    p.add_level("redefine", process_redefine)

    def process_retype(tokens):
        nonlocal last_reaction
        n_reac = len(last_reaction.reactants)
        if n_reac == 0:
            raise ParseError("[reactants] must come before [retype]")

        i, j = parse_pair(tokens, 0, "pair")
        ntype = unwrap(tokens, 1, "word")
        last_reaction.retypes.append((i, j, ntype))

    p.add_level("retype", process_retype)

    def process_rename(tokens):
        nonlocal last_reaction
        n_reac = len(last_reaction.reactants)
        if n_reac == 0:
            raise ParseError("[reactants] must come before [rename]")

        i, j = parse_pair(tokens, 0, "pair")
        nname = unwrap(tokens, 1, "word")
        last_reaction.renames.append((i, j, nname))

    p.add_level("rename", process_rename)

    def process_remass(tokens):
        nonlocal last_reaction
        n_reac = len(last_reaction.reactants)
        if n_reac == 0:
            raise ParseError("[reactants] must come before [remass]")

        i, j = parse_pair(tokens, 0, "pair")
        nmass = unwrap(tokens, 1, "float")
        last_reaction.remasses.append((i, j, nmass))

    p.add_level("remass", process_remass)

    def process_recharge(tokens):
        nonlocal last_reaction
        n_reac = len(last_reaction.reactants)
        if n_reac == 0:
            raise ParseError("[reactants] must come before [recharge]")

        i, j = parse_pair(tokens, 0, "pair")
        ncharge = unwrap(tokens, 1, "float")
        last_reaction.recharges.append((i, j, ncharge))

    p.add_level("recharge", process_recharge)

    def process_soft_core(tokens):
        nonlocal last_reaction
        n_reac = len(last_reaction.reactants)
        if n_reac == 0:
            raise ParseError("[reactants] must come before [soft_core]")

        i, j = parse_pair(tokens, 0, "pair")
        sc_lam = unwrap(tokens, 1, "float")
        sc_alpha = unwrap(tokens, 2, "float")
        last_reaction.soft_core.append((i, j, sc_lam, sc_alpha))

    p.add_level("soft_core", process_soft_core)

    def process_system(tokens):
        nonlocal system_defined
        require_complete_reaction()
        system_defined = True  # for a hack so that reactions are complete

    p.add_level("system", process_system, mandatory=True, unique=True)

    # run parser
    ok = p.parse(file, include_dir, defines)
    if not ok:
        return False, None
    else:
        return True, (system, topology)
