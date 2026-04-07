import numpy as np
from ..__rust import Fragment, FragList, DetectionTemplate, DetectionTemplateList, detection, PeriodicBox
from .modification_template import ModificationTemplate
from .graph import Graph, GraphMatch, match_atoms
from ..__core import System

class TopStar:
    """
    The glue between the following components:
    - graph matching algorithm
    - detection algorithm
    - modification algorithm
    """
    def __init__(self, system: System):
        self.system = system
        self.n_atoms = system.atom_count()
        self.frag_list = FragList(self.n_atoms)
        self.detection_templates = DetectionTemplateList()
        for d in system.additional_data.get("detection_templates", []):
            self.detection_templates.add_detection_template(d)
        self.graphs: dict[str, Graph] = system.additional_data.get("graphs") or {}

        # run the graph matching algorithm on all molecules separately
        i = 0
        for mol, n in system.initial_molecules:
            n_atoms = len(system.molecule_types.get(mol).atoms)
            for j in range(n):
                self.try_match_graphs(
                    set(range(i, i + n_atoms))
                )
                i += n_atoms

    def try_match_graphs(self, atoms: set[int]) -> None:
        """
        Given a set of atoms, find all graph matches of all known
        graphs and add them to the fragment list.

        Should be called after system is constructed.

        Initial system construction: call this for every molecule.

        Reaction: Atoms should be all affected atoms in a reaction,
        and all their neighbors. Remove all pre-existing fragments
        that contain atoms on which try_match_graphs is called first!
        """

        matches: list[GraphMatch] = []
        for graph in self.graphs.values():
            matches += match_atoms(
                graph, atoms, self.system
            )
        for m in matches:
            frag_atoms = []
            for key, _, _, _ in m.graph.atoms:
                val = m.atoms.get(key)
                # key - name in the graph
                # val - atom id
                if val is None:
                    frag_atoms.append(-1)
                else:
                    frag_atoms.append(val)

            self.frag_list.add_fragment(
                m.graph.name, frag_atoms
            )

    def update_rates(self, reactions: list[tuple[str, list[int]]], counts: dict[str, int], volume: float) -> None:
        """
        Updates "observed rate" in detection templates based on the reactions happening and current
        reactant concentrations (specified using counts and volume).
        """
        # TODO
        # note we only have post-rate adjusted rates
        pass

    def detection(self, pbc: PeriodicBox, pos) -> list[tuple[str, list[int]]]:
        """
        Runs the detection algorithm, given a periodic box and atom positions and current state in TopStar.
        Returns the list of reactions.
        """
        return detection(
            self.frag_list,
            self.detection_templates,
            pbc,
            pos
        )

    def modification(self, reactions: list[tuple[str, list[int]]]) -> list[tuple[str, list[Fragment]]]:
        """
        Runs the modification algorithm.

        Modifies TopStar and System according to the reaction templates.

        Returns the list of reactions that were successfully applied to the system.
        """
        completed = []

        for (rx, frag_ids) in reactions:
            frags = [
                self.frag_list.get_fragment(frag_id)
                for frag_id in frag_ids
            ]
            if any(frag is None for frag in frags):
                # pass reactions if a previous reactions' modification algorithm destroyed the reactant fragment
                # of another reaction
                continue

            flattened_atoms = []
            for f in frags:
                flattened_atoms.extend(f.atoms)

            # execute modification tempate
            m_template = self.system.molecule_types[rx]
            assert isinstance(m_template, ModificationTemplate)
            m_template.instantiate(self.system, flattened_atoms)

            # remove old graphs
            recalc = self.system.populate_neighbors(
                x for x in flattened_atoms if x >= 0
            )
            self.frag_list.delete_fragments_for_atoms(list(recalc))

            # add new graphs
            self.try_match_graphs(recalc)

            # mark reaction as completed
            # return frags, since frag_id's now refer to non existent frags
            completed.append((rx, frags))

        return completed

    def toggle_softcore(self, reactions: list[tuple[str, list[Fragment]]], on: bool):
        for (rx, frags) in reactions:
            m_template = self.system.molecule_types[rx]
            flattened_atoms = []
            for f in frags:
                flattened_atoms.extend(f.atoms)

            m_template.toggle_softcore(self.system, flattened_atoms, on)

