from __future__ import annotations

import math
from typing import TextIO

import numpy as np
import openmm as mm

from ..__rust import Fragment
from ..__simulation import Reporter, Simulation


class GlobalMinimizer(Reporter):
    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        pass

    def on_simulation_finish(self, simulation: Simulation) -> None:
        pass

    def __init__(
        self,
        minimization_steps: int = 500,
        _restraint_force_global: float = 1000,
        _restraint_force_local: float = 50,
        r_movable: float = 1.0,
        whole_molecule: bool = True,
        harmonic_constraints: bool = True,
    ) -> None:
        """
        Global Minimizer. Uses the Reporter API to minimize the system after a reaction.

        :param minimization_steps: The number of steps to run the minimization algorithm for.
        :param restraint_force_global: Strength of position restraints for the entire system. (TODO)
        :param restraint_force_local: Strength of position restraint for reacting atoms and their
            direct neighbors. (TODO)
        :param r_movable: Distance cutoff of how far from reacting atoms counts as local, in nanometers. Set to 0. to only consider the bond graph.
        :param whole_molecule: If True, the entire molecule will be considered local and included in the minimization. If false, only reactant atoms and their bond neighbors will be included.
        :param harmonic_constraints: If True, the constraints will be converted to stiff harmonic bonds during minimization.
        """
        assert minimization_steps > 0, "Must specify minimization_steps > 0"
        self.minimization_steps = minimization_steps
        self.r_movable = r_movable
        self.whole_molecule = whole_molecule
        self.harmonic_constraints = harmonic_constraints

    def pre_simulation_start(self, simulation: Simulation) -> None:
        assert simulation.integrator is not None

    def on_reaction(
        self, simulation: Simulation, reactions: list[tuple[str, list[Fragment]]]
    ) -> None:
        # save velocities
        simulation.info("on_reaction Global Minimizer")

        simulation.top.toggle_softcore(reactions, True)
        if self.harmonic_constraints:
            simulation.system.toggle_constraints_as_harmonic_bonds(True)

        pos, box = simulation.context.get_positions()
        # set movable
        # atoms in reacting atoms only
        atoms = set()
        for _, frags in reactions:
            for frag in frags:
                for atom in frag.atoms:
                    if atom != -1:
                        atoms.add(atom)

        if self.r_movable > 0.0:
            atoms |= box.which_atoms_within_distance(pos, atoms, self.r_movable)

        atoms = simulation.system.populate_neighbors(
            atoms, recursive=self.whole_molecule
        )
        simulation.info("atoms chosen")

        # TODO constraints

        simulation.context.minimize_energy(max_steps=self.minimization_steps)

        if self.harmonic_constraints:
            simulation.system.toggle_constraints_as_harmonic_bonds(False)
        simulation.top.toggle_softcore(reactions, False)
        simulation.info("preparing to set integrator back")
