def _():
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
