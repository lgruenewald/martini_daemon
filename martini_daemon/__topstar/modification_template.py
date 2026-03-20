from ..__core import MoleculeType
from ..__parser import TokenList, ParseException, TokenParseException, GromacsTopFile


class ModificationTemplate(MoleculeType):
    def __init__(self, parent: GromacsTopFile, path: str, line_num: int):
        super().__init__(parent, path, line_num)
        self.topology = parent.result.get("topology")
        self.__reactants = []
        # TODO ...

    def add_reactant(self):
        # TODO verify that reactant exists before adding it to reactants[]
        pass

    def parse_index(self, tokens: TokenList, index: int) -> int:
        n_reactants = len(self.__reactants)
        reactant_index, atom_name = tokens.unwrap(index, "pair")
        # == Verification ==
        if reactant_index < 0 or reactant_index >= n_reactants:
            raise TokenParseException(
                tokens[index],
                f"Reactant index {reactant_index+1} out of range: 1 to {n_reactants}."
            )
        graph_name = self.__reactants[reactant_index]
        # verified on add_reactant
        graph = self.topology.graphs[graph_name]
        atom_index = graph.atom_name_to_index.get(atom_name)
        if atom_index is None:
            raise TokenParseException(
                tokens[index],
                f"Atom {atom_name} not found in reactant {graph_name}."
            )
        if graph.atoms[atom_index][3] not in {
            GraphAtomType.NORMAL, GraphAtomType.OPT
        }:
            raise TokenParseException(
                tokens[index],
                "Attempt to reference a forbidden atom!"
            )
        # get index
        flattened_index = sum(
            len(self.topology.graphs[self.__reactants[i]].atoms)
            for i in range(reactant_index)
        ) + atom_index
        return flattened_index

    def add_atoms_to_system(self, system):
        # reactions don't add atoms to system any more
        raise ParseException(
            "Attempt to add a reaction to a [system]/[molecules]. Please use a [moleculetype] molecule."
        )
