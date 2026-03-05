from .gromacs_top_format import directive

@directive
class Molecule:
    pass
# TODO add a "process_nrel_excl" here instead of adding exclusions together with bonds