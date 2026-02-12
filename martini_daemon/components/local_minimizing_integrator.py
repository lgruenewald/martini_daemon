import openmm as mm
import numpy as np
from .daemon_integrator import DaemonIntegrator
import math


class LocalMinimizingIntegrator(DaemonIntegrator):
    """
        Inspired by
        https://github.com/choderalab/openmmtools/blob/main/openmmtools/integrators.py
        GradientDescentMinimizationIntegrator
    """

    def __init__(
        self, dt_ps, T_K, friction_ps1, minimizer, minimization_steps=500,
        report_every=0, check_convergence_every=10
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
        self.check_convergence_every = check_convergence_every

    def reset(self, shape):
        self.minimizer.reset(shape)

    def report(self, i, rem, op="a"):
        with open(f"{self.sim_name}_minimization{i}.log", op) as f:
            f.write(f"Simulation frame {i}, remaining steps {rem}\n")
            f.write("==============================================\n")
            f.write(self.minimizer.report())

    def set_reactions(self, reactions, system, top, i):
        # save vels
        vels = system._context.getState(
            velocities=True
        ).getVelocities(
            asNumpy=True
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

        if self.report_every > 0:
            gcd = math.gcd(self.report_every, self.check_convergence_every)
        else:
            gcd = self.check_convergence_every
        remaining = self.minsteps

        if self.report_every > 0:
            self.report(i, remaining, op="w")

        while remaining > 0:
            csteps = min(gcd, remaining)
            self.integrator.step(csteps)
            remaining -= csteps
            if self.report_every > 0 and remaining % self.report_every == 0:
                self.report(i, remaining)
            # check convergence every is guaranteed to check AT LEAST as often as it says
            if self.minimizer.getGlobalVariableByName("converged") == 1:
                break
            
        # reporters and cleanup
        # TODO find out why velocities change significantly during minimization
        system.set_velocities(vels)
        for rep in self.reporters:
            rep.post_di_minimize()
        self.integrator.setCurrentIntegrator(0)

    def step(self, n_steps):
        self.integrator.step(n_steps)

    def get_integrator(self):
        return self.integrator

    def getStepSize(self):
        return self.dt * mm.unit.picosecond
