from __future__ import annotations

import openmm as mm

from ..__core import Force
from ..__rust import Fragment
from ..__simulation import Reporter, Simulation


class MinimizationRestraint(Force):
    def delta_degrees_of_freedom(self) -> int:
        return 0

    @classmethod
    def is_coupling(cls) -> bool:
        return False

    @classmethod
    def get_name(cls) -> str:
        return "minimization_restraint"

    def _set_force_obj(self) -> mm.Force:
        # FIXME: verify non 90 degree angles
        force = mm.CustomExternalForce(
            "0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2;"
        )
        force.addPerParticleParameter("k")
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")

        for i in range(self.system.num_atoms()):
            force.addParticle(i, [0., 0., 0., 0.])

        return force

    def __update_context(self) -> None:
        """Update context without reinitializing (FIXME: check why this is slow)."""
        if self._force is not None:
            assert isinstance(self._force, mm.CustomExternalForce)
            if (ctx := self.system._get_context()) is not None:
                self._force.updateParametersInContext(ctx)

    def reset(self) -> None:
        """Set all restraints to 0, update context without reinitialize."""
        if self._force is not None:
            assert isinstance(self._force, mm.CustomExternalForce)
            for i in range(self.system.num_atoms()):
                self._force.setParticleParameters(i, i, [0., 0., 0., 0.])
            self.system.flag_reinitialize()
        # FIXME switch to self.update_context()

    def set_atom_restraints(self, i: int, k: float, x: float, y: float, z: float) -> None:
        """
        Set parameters for a single atom.

        Note: must only call if the force has been built.
        """
        assert self._force is not None
        assert isinstance(self._force, mm.CustomExternalForce)
        self._force.setParticleParameters(i, i, [k, x, y, z])
        self.system.flag_reinitialize()
        # FIXME rely on self.__update_context()

class GlobalMinimizer(Reporter):
    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        pass

    def on_simulation_finish(self, simulation: Simulation) -> None:
        pass

    def __init__(
        self,
        minimization_steps: int = 50,
        restraint_force_global: float = 1000,
        restraint_force_local: float = 50,
        r_movable: float = 0.0,
        whole_molecule: bool = False,
    ) -> None:
        """
        Global Minimizer. Uses the Reporter API to minimize the system after a reaction.

        :param minimization_steps: The number of steps to run the minimization algorithm for.
        :param restraint_force_global: Strength of position restraints for the entire system.
        :param restraint_force_local: Strength of position restraint for reacting atoms and their
            direct neighbors.
        :param r_movable: Distance cutoff of how far from reacting atoms counts as local, in nanometers. Set to 0. to only consider the bond graph.
        :param whole_molecule: If True, the entire molecule will be considered local and included in the minimization. If false, only reactant atoms and their bond neighbors will be included.
        """
        assert minimization_steps > 0, "Must specify minimization_steps > 0"
        self.minimization_steps = minimization_steps
        self.r_movable = r_movable
        self.whole_molecule = whole_molecule
        self.restraint_force = None
        self.k_local = restraint_force_local
        self.k_global =restraint_force_global

    def pre_simulation_start(self, simulation: Simulation) -> None:
        assert simulation.integrator is not None
        self.restraint_force = MinimizationRestraint(simulation.system)
        simulation.system.add_force(self.restraint_force)

    def on_reaction(
        self, simulation: Simulation, reactions: list[tuple[str, list[Fragment]]]
    ) -> None:
        # save velocities
        simulation.info("on_reaction Global Minimizer")
        assert self.restraint_force is not None

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

        for i in range(simulation.system.num_atoms()):
            k = self.k_local if i in atoms else self.k_global
            self.restraint_force.set_atom_restraints(i, k, pos[i, 0], pos[i, 1], pos[i, 2])

        simulation.info("restraints set")
        simulation.info("minimization start")
        # minimization should reinitialize anyway -- there were reactions!
        simulation.context.minimize_energy(max_steps=self.minimization_steps)

        simulation.info("minimization finished")
        self.restraint_force.reset()
        simulation.info("restraints reset")
