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
        report_every=0
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
        self.report_every = report_every

    def reset(self, shape):
        self.minimizer.reset(shape)

    def report(self, i, rem):
        with open(f"{self.sim_name}_minimization{i}.log", "a") as f:
            f.write(f"Simulation frame {i}, remaining steps {rem}\n")
            f.write("==============================================\n")
            f.write(self.minimizer.report())

    def set_reactions(self, reactions, system, top, i):
        # save vels and zero them out
        vels = system._context.getState(
            velocities=True
        ).getVelocities(
            asNumpy=True
        )
        system._context.setVelocities(
            np.zeros(shape=vels.shape)
        )
        self.integrator.setCurrentIntegrator(1)
        # integrator state setup
        self.reset(vels.shape)
        # set movable
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
        if self.report_every == 0:
            self.integrator.step(self.minsteps)
        else:
            # create/empty file
            with open(f"{self.sim_name}_minimization{i}.log", "w") as f:
                f.write("")
            self.report(i, self.minsteps)
            remaining = self.minsteps
            while remaining > 0:
                csteps = min(self.report_every, remaining)
                self.integrator.step(csteps)
                remaining -= csteps
                self.report(i, remaining)
        # reporters and cleanup
        system._context.setVelocities(vels)
        for rep in self.reporters:
            rep.post_di_minimize()
        self.integrator.setCurrentIntegrator(0)

    def step(self, n_steps):
        self.integrator.step(n_steps)

    def get_integrator(self):
        return self.integrator

    def getStepSize(self):
        return self.dt * mm.unit.picosecond
