
class DaemonIntegrator():
    reporters = None

    def set_reactions(self, reactions, system, top, i):
        pass

    def step(self, n_steps):
        pass

    def get_integrator(self):
        pass

    def getStepSize(self):
        pass

    def add_reporter(self, rep):
        if self.reporters is None:
            self.reporters = [rep]
        else:
            self.reporters.append(rep)
