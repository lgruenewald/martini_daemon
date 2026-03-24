def _():

    def instantiate(self, molname: str) -> None:
        """
        Adds a fresh new copy of moleculetype molname to the system.
        R
        """
        # molname must refer to a Molecule type
        mol = self.molecules.get(molname)
        if mol is None:
            raise ValueError(f"Can't find mol {molname}")

        # add atoms to S* and other bookkeeping
        atoms = []
        prev_resnum = 0
        for atom in mol.atoms:
            type, resnum, resname, atomname, chargegr, charge, mass = atom
            if resnum != prev_resnum:
                self.system.new_residue()
                prev_resnum = resnum
            p = self.system.add_atom(atomname, resname, type, charge, mass, 1, 0.5)
            atoms.append(p)
            self.defrag_list.append([])
            self.interaction_list.append([])

        # instantiate interactions
        self.instantiate_over_existing(mol, atoms)

        # add graphs to system
        self.try_match_graphs(set(atoms), molname)

        # do initial molecules info, used e.g. in old_helpers/monomer
        if (
            len(self.initial_molecules) > 0
            and self.initial_molecules[-1][0] == molname
        ):
            _, n, n_atoms = self.initial_molecules[-1]
            assert len(atoms) == n_atoms
            self.initial_molecules[-1] = (molname, n+1, n_atoms)
        else:
            self.initial_molecules.append((molname, 1, len(atoms)))

    def instantiate_over_existing(
        self, mol: Molecule, atoms: list[Fragment] | list[int],
    ) -> None:
        """
        Takes a name of a molecule, adds interactions to those atoms
        according to the molecule.
        """
        def index_pair(index: int | tuple[int, int]) -> int:
            # abstraction, because of the duality of Molecule
            if type(index) is int:
                return atoms[index]
            else:
                idi, atomi = index
                return atoms[idi].atoms[atomi]

        # exclusions
        for (i, j) in mol.exclusions:
            pi = index_pair(i)
            pj = index_pair(j)
            if pi == -1 or pj == -1:
                # optional atom missing
                continue
            if pi < pj:
                e = self.system.exclusions.add(pi, pj)
                self.interaction_list[pi].append(e)
                self.interaction_list[pj].append(e)
        # posres
        for (i, kx, ky, kz) in mol.posres:
            if self.respos is None:
                raise ValueError(
                    "respos is None but there are position restraints."
                )
            pi = index_pair(i)
            x0, y0, z0 = self.respos[pi]
            self.system.posres.add((pi), (kx, ky, kz, x0, y0, z0))
        # generic interactions
        for (force, members, params) in mol.interactions:
            member_atoms = []
            missing_opt = False
            for x in members:
                atom = index_pair(x)
                if atom == -1:
                    missing_opt = True
                    break
                member_atoms.append(atom)
            if missing_opt:
                # optional missing => skip
                continue
            f = force.add(member_atoms, params)
            for member in member_atoms:
                self.interaction_list[member].append(f)

    def update_frag_counts(self):
        """
        Sets self.frag_counts. Call after every modification.

        self.frag_counts is used notably for rate calculation.
        """
        # TODO too much state sync..
        self.frag_counts = {}
        for v in self.frag_list.values():
            k = v.name
            if self.frag_counts.get(k) is None:
                self.frag_counts[k] = 1
            else:
                self.frag_counts[k] += 1

    # ======= (2/3) Detection things - see old_detection.pyx =======
    def init_dm(self, name) -> None:
        # TODO remove this function, reduce fragility that way
        # called exactly once after parsing or loading from file is finished
        self.update_frag_counts()
        self.out_name = name
        for reporter in self.reporters:
            reporter.init_dm(name)

    def update_observed_rate(self, rx: ReactionTemplate):
        # update_observed_rate only called if it is defined => not none
        # relative_rate is parsed as "positive" => non zero, safe to divide
        rate = rx.reaction_counter / rx.relative_rate
        for reactant in rx.reactants:
            if self.frag_counts.get(reactant) in {0, None}:
                # no reactant => no reaction
                # keeps old observed_rate and we avoid divisions by 0
                return
            rate /= self.frag_counts[reactant]

        # if it's not initialized yet make an initial value
        if rx.observed_rate is None:
            rx.observed_rate = rate
        else:
            # otherwise smooth it
            rx.observed_rate = rate * 0.01 + rx.observed_rate * 0.99

        predicted = rx.observed_rate * self.highest_probability

        # if the prediction value goes below 0, 0 it and issue a warning
        if predicted < 0.:
            self.logger.warn(
                f"The smoothed rate for reaction {rx.name} went below 0."
            )
            rx.observed_rate = 0.
            predicted = 0.

        # set the relative rate = 1 value to the slowest reaction
        if self.absolute_rate is None or predicted < self.absolute_rate:
            self.absolute_rate = predicted

    def detection(
        self, step: int, box, pos
    ) -> list[tuple[list[Fragment], ReactionTemplate]]:
        for reporter in self.reporters:
            reporter.pre_detection(step, self.out_name)
        rxs = []
        for more_rxs in self.reactions.values():
            for rx in more_rxs:
                rxs.append(rx)
        reactions = detection_all(
            self, box, pos
        )
        """
        reactions = detection_all(
            self.frag_list, rxs,
            box, pos,
            self.nlist_cutoff,
            self.absolute_rate,
            # TODO save state (not sure there is much point)
            np.random.default_rng()
        )
        """
        # preparations for the next step
        self.absolute_rate = self.max_absolute_rate
        for _, rxs in self.reactions.items():
            for rx in rxs:
                if rx.relative_rate is not None:
                    self.update_observed_rate(rx)
                rx.reaction_counter = 0
        if self.absolute_rate is None:
            self.absolute_rate = -1.
        return reactions

    # ======= (3/3) Modification things =======
    def clean_defrag(self, frag: Fragment, atom: int) -> None:
        """Removes references to frag in the defrag_list for atom_id atom."""
        self.defrag_list[atom] = list(filter(
            lambda x: x != frag.frag_id,
            self.defrag_list[atom]
        ))

    def remove_fragment(self, frag: Fragment) -> None:
        """Removes a fragment from frag_list and defrag_list"""

        for atom in frag.atoms:
            if atom != -1:
                self.clean_defrag(frag, atom)

        del self.frag_list[frag.frag_id]

    def remove_interaction(self, interaction: Interaction) -> None:
        """Removes an interaction from interaction_list and S*"""
        for atom in interaction.get_members():
            self.interaction_list[atom] = list(filter(
                lambda x: x != interaction,
                self.interaction_list[atom]
            ))
        interaction.remove()

    def process_break(self, frags: list[Fragment], rx: ReactionTemplate,
                      ) -> None:
        """Process [break] in rx over frags."""

        for group in rx.break_groups:
            group_atoms = []
            all_found = True
            for id, atom in group:
                atom = frags[id].atoms[atom]
                if atom == -1:
                    all_found = False
                    break
                group_atoms.append(atom)
            if not all_found:
                # missing optional
                continue
            # breaking interactions
            # for every group in break_groups
            # check all interactions and break if all in group are in inter
            for inter in self.interaction_list[group_atoms[0]][:]:
                if all(map(
                           lambda group_member:
                           group_member in inter.get_members(),
                           group_atoms
                       )):
                    self.remove_interaction(inter)

    def process_update(self, frags: list[Fragment], rx: ReactionTemplate,
                       ) -> None:
        """Process [rx_update] in rx over frags."""

        for group in rx.update_groups:
            group_atoms = []
            all_found = True
            for id, atom in group:
                atom = frags[id].atoms[atom]
                if atom == -1:
                    all_found = False
                    break
                group_atoms.append(atom)
            if not all_found:
                # missing optional
                continue
            # breaking interactions
            # check that all in interaction are in group
            for group_atom in group_atoms:
                for inter in self.interaction_list[group_atom][:]:
                    if all(map(
                               lambda inter_member: inter_member in group_atoms,
                               inter.get_members()
                           )):
                        self.remove_interaction(inter)

    def populate_neighbors(self, atoms: set[int]) -> set[int]:
        """
        For a set of atoms, return a set that also contains their neighbors.

        Neighbor = shared interaction (e.g. bond, angle, exclusion...)
        """
        res = set()
        for atom in atoms:
            res.add(atom)
            for inter in self.interaction_list[atom]:
                for member in inter.get_members():
                    res.add(member)
        return res

    def remove_overlapping_graphs(self, atoms: set[int]) -> None:
        """Remove all fragments in T* that contain any of the given atoms."""
        for atom in atoms:
            for frag_id in self.defrag_list[atom][:]:
                frag = self.frag_list.get(frag_id)
                if frag is not None and frag.graph is not None:
                    self.remove_fragment(frag)

    def template_update_atoms(
        self, frags: list[Fragment], rx: ReactionTemplate
    ) -> None:
        """Process [rename], [retype], [recharge], [remass]."""
        # re* overrides everything
        for (id, atom, new_name) in rx.renames:
            atom_id = frags[id].atoms[atom]
            if atom_id == -1:
                continue
            self.system.rename(atom_id, new_name)
        for (id, atom, new_type) in rx.retypes:
            atom_id = frags[id].atoms[atom]
            if atom_id == -1:
                continue
            self.system.retype(atom_id, new_type)
        for (id, atom, new_charge) in rx.recharges:
            atom_id = frags[id].atoms[atom]
            if atom_id == -1:
                continue
            self.system.recharge(atom_id, new_charge)
        for (id, atom, new_mass) in rx.remasses:
            atom_id = frags[id].atoms[atom]
            if atom_id == -1:
                continue
            self.system.remass(atom_id, new_mass)

    def modification(
        self,
        i: int,
        reactions: list[tuple[list[Fragment], ReactionTemplate]]
    ) -> None:
        """Modification helper for the D/M algorithm
        """

        for reporter in self.reporters:
            reporter.pre_modification(i, reactions, self.out_name)

        completed_reactions = []

        for (frags, rx) in reactions:
            # check if all reactants still exist
            if any(
                [self.frag_list.get(frag.frag_id) is None for frag in frags]
            ):
                # we have no way of detecting
                # this in the detection algorithm, since this is
                # about graphs being recalculated during reactions
                # potentially invalidating other graphs that react
                # in the same frame
                continue

            # just normal atoms to instantiate products over
            product_atoms = []
            # which atoms to recalculate graphs over
            graph_recalc = set()
            for f in frags:
                product_atoms += list(filter(lambda x: x != -1, f.atoms))
                graph_recalc |= set(filter(lambda x: x != -1, f.atoms))
            # add neighbors since those can be changed too (opt/not atoms)
            graph_recalc = self.populate_neighbors(set(product_atoms))

            # [rx_break]
            self.process_break(frags, rx)
            # [rx_update]
            self.process_update(frags, rx)

            # graphs get recalculated later over the same atoms
            self.remove_overlapping_graphs(graph_recalc)

            self.template_update_atoms(frags, rx)
            self.instantiate_over_existing(self.molecules[rx.name], frags)

            self.try_match_graphs(graph_recalc)

            completed_reactions.append((frags, rx))

        self.update_frag_counts()

        for reporter in self.reporters:
            reporter.post_modification(i, self.out_name, completed_reactions)

        return completed_reactions

    def toggle_sc(
            self, reactions: list[tuple[list[Fragment], ReactionTemplate]],
            toggle: bool
    ) -> None:
        for (frags, rx) in reactions:
            for (id, atom, new_lam, new_alpha) in rx.soft_core:
                atom_id = frags[id].atoms[atom]
                if atom_id == -1:
                    continue
                if toggle:
                    self.system.update_sc(atom_id, new_lam, new_alpha)
                else:
                    self.system.update_sc(atom_id, 1, 0.5)

#            for f in frags:
#                for atom_id in f.atoms:
#                    if atom_id == -1:
#                         continue
#                    if toggle:
#                        self.system.update_sc(atom_id, 0.6, 0.5)
#                    else:
#                        self.system.update_sc(atom_id, 1, 0.5)
        return

