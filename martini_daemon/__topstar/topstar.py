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
        self.absolute_rate = None

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

    def detection(self, pbc: PeriodicBox, pos) -> list[tuple[str, list[int]]]:
        return detection(
            self.frag_list,
            self.detection_templates,
            self.absolute_rate,
            pbc,
            pos
        )