
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
        if len(last_reaction.__reactants) > 0:
            raise ParseError("Only one set of reactants per reaction.")
        last_reaction.__reactants = rxs

    p.add_level("reactants", process_reactants)
    # alias:
    p.add_level("reactant", process_reactants)

    def process_conditions(tokens):
        nonlocal last_reaction
        n_reac = len(last_reaction.__reactants)
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
            case "probability":
                last_reaction.probability = unwrap(tokens, 1, "positive")
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
        n_reac = len(last_reaction.__reactants)
        if n_reac == 0:
            raise ParseError("[reactants] must come before [break]")

        atoms = []
        for i in range(len(tokens)):
            atoms.append(parse_pair(tokens, i, "pair"))
        last_reaction.break_groups.append(atoms)

    p.add_level("break", process_break)

    def process_update(tokens):
        nonlocal last_reaction
        n_reac = len(last_reaction.__reactants)
        if n_reac == 0:
            raise ParseError("[reactants] must come before [update]")

        atoms = []
        for i in range(len(tokens)):
            atoms.append(parse_pair(tokens, i, "pair"))
        last_reaction.update_groups.append(atoms)

    p.add_level("update", process_update)

    def process_redefine(tokens):
        nonlocal last_reaction
        n_reac = len(last_reaction.__reactants)
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
        n_reac = len(last_reaction.__reactants)
        if n_reac == 0:
            raise ParseError("[reactants] must come before [retype]")

        i, j = parse_pair(tokens, 0, "pair")
        ntype = unwrap(tokens, 1, "word")
        last_reaction.retypes.append((i, j, ntype))

    p.add_level("retype", process_retype)

    def process_rename(tokens):
        nonlocal last_reaction
        n_reac = len(last_reaction.__reactants)
        if n_reac == 0:
            raise ParseError("[reactants] must come before [rename]")

        i, j = parse_pair(tokens, 0, "pair")
        nname = unwrap(tokens, 1, "word")
        last_reaction.renames.append((i, j, nname))

    p.add_level("rename", process_rename)

    def process_remass(tokens):
        nonlocal last_reaction
        n_reac = len(last_reaction.__reactants)
        if n_reac == 0:
            raise ParseError("[reactants] must come before [remass]")

        i, j = parse_pair(tokens, 0, "pair")
        nmass = unwrap(tokens, 1, "float")
        last_reaction.remasses.append((i, j, nmass))

    p.add_level("remass", process_remass)

    def process_recharge(tokens):
        nonlocal last_reaction
        n_reac = len(last_reaction.__reactants)
        if n_reac == 0:
            raise ParseError("[reactants] must come before [recharge]")

        i, j = parse_pair(tokens, 0, "pair")
        ncharge = unwrap(tokens, 1, "float")
        last_reaction.recharges.append((i, j, ncharge))

    p.add_level("recharge", process_recharge)

    def process_soft_core(tokens):
        nonlocal last_reaction
        n_reac = len(last_reaction.__reactants)
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
