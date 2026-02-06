import openmm as mm
from .daemon_integrator import DaemonIntegrator

class ReactionSensitiveLangevinIntegrator(DaemonIntegrator):

    def __init__(
        self, dt_ps, T_K, friction_ps1, subdivision=4, equilibration_length=10,
        minimization_steps=0
    ):
        self.dt = dt_ps
        self.T = T_K
        self.friction = friction_ps1
        self.subdivision = subdivision
        self.eqlen = equilibration_length
        self.integrator = mm.CompoundIntegrator()
        self.integrator.addIntegrator(
            mm.LangevinMiddleIntegrator(
                T_K, friction_ps1, dt_ps
            )
        )
        self.integrator.addIntegrator(
            mm.LangevinIntegrator(
                T_K, friction_ps1, dt_ps / subdivision
            )
        )
        self.remaining = 0
        self.minsteps = minimization_steps

    def set_reactions(self, reactions, system, top):
        self.remaining = self.eqlen
        self.integrator.setCurrentIntegrator(1)
        if self.minsteps > 0:
            state = system._context.getState(
                velocities=True
            )
            system.minimize_energy(
                max_steps=self.minsteps
            )
            system._context.setVelocities(
                state.getVelocities(
                    asNumpy=True
                )
            )
            if self.reporters is not None:
                for rep in self.reporters:
                    rep.post_di_minimize()

    def finish_equilibration(self):
        self.integrator.setCurrentIntegrator(0)
        if self.reporters is not None:
            for rep in self.reporters:
                rep.post_di_equilibrate()

    def step(self, n_steps):
        if self.remaining >= n_steps:
            self.remaining -= n_steps
            self.integrator.step(n_steps * self.subdivision)
            if self.remaining == 0:
                self.finish_equilibration()
        elif self.remaining > 0:
            self.integrator.step(self.remaining * self.subdivision)
            self.finish_equilibration()
            self.integrator.step(n_steps - self.remaining)
            self.remaining = 0
        else:
            self.integrator.step(n_steps)

    def get_integrator(self):
        return self.integrator

    def getStepSize(self):
        return self.dt * mm.unit.picosecond
