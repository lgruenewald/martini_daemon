def _():

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
