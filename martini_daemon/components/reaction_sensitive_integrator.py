import openmm as mm


class DaemonIntegrator():
    def set_reactions(self, reactions):
        pass

    def step(self, n_steps):
        pass

    def get_integrator(self):
        pass

    def getStepSize(self):
        pass


class ReactionSensitiveLangevinIntegrator(DaemonIntegrator):

    def __init__(
        self, dt_ps, T_K, friction_ps1, subdivision=4, equilibration_length=50
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

    def set_reactions(self, reactions):
        self.remaining = self.eqlen
        self.integrator.setCurrentIntegrator(1)

    def step(self, n_steps):
        if self.remaining >= n_steps:
            self.remaining -= n_steps
            self.integrator.step(n_steps * self.subdivision)
            if self.remaining == 0:
                self.integrator.setCurrentIntegrator(0)
        elif self.remaining > 0:
            self.integrator.step(self.remaining * self.subdivision)
            self.integrator.setCurrentIntegrator(0)
            self.integrator.step(n_steps - self.remaining)
            self.remaining = 0
        else:
            self.integrator.step(n_steps)

    def get_integrator(self):
        return self.integrator

    def getStepSize(self):
        return self.dt * mm.unit.picosecond
