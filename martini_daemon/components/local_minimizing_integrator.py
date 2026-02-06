import openmm as mm
import numpy as np
from .daemon_integrator import DaemonIntegrator


class LocalMinimizingIntegrator(DaemonIntegrator):
    """
        Inspired by
        https://github.com/choderalab/openmmtools/blob/main/openmmtools/integrators.py
        GradientDescentMinimizationIntegrator
    """

    def __init__(
        self, dt_ps, T_K, friction_ps1, minimizer, minimization_steps=500,
    ):
        self.dt = dt_ps
        self.T = T_K
        self.friction = friction_ps1
        self.integrator = mm.CompoundIntegrator()
        self.integrator.addIntegrator(
            mm.LangevinMiddleIntegrator(
                T_K, friction_ps1, dt_ps
            )
        )
        self.minimizer = minimizer
        self.integrator.addIntegrator(minimizer)
        assert minimization_steps > 0, "Must specify minimization_steps > 0"
        self.minsteps = minimization_steps

    def set_reactions(self, reactions, system, top):
        # save vels and zero them out
        vels = system._context.getState(
            velocities=True
        ).getVelocities(
            asNumpy=True
        )
        system._context.setVelocities(
            np.zeros(shape=vels.shape)
        )
        # set movable
        self.integrator.setCurrentIntegrator(1)
        movable = np.zeros(shape=vels.shape)
        atoms = set()
        for (frags, _) in reactions:
            for frag in frags:
                for atom in frag.atoms:
                    if atom != -1:
                        atoms.add(atom)
        for atom in atoms.copy():
            for inter in top.interaction_list[atom]:
                for member in inter.get_members():
                    atoms.add(member)
        for atom in atoms:
            movable[atom, :] = 1.
        self.minimizer.setPerDofVariableByName(
            "movable",
            movable
        )
        # minimize
        self.integrator.step(self.minsteps)
        # reporters and cleanup
        system._context.setVelocities(vels)
        for rep in self.reporters:
            rep.post_di_minimize()
            rep.post_di_equilibrate()
        self.integrator.setCurrentIntegrator(0)

    def step(self, n_steps):
        self.integrator.step(n_steps)

    def get_integrator(self):
        return self.integrator

    def getStepSize(self):
        return self.dt * mm.unit.picosecond
