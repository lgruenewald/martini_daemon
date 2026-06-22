import openmm as mm

from ..__rust import Fragment
from ..__simulation import Reporter, Simulation


class GlobalIntegratorMinimizer(Reporter):
    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        pass

    def on_simulation_finish(self, simulation: Simulation) -> None:
        pass

    def __init__(
        self,
        minimizer: mm.Integrator,
        n_steps: int = 10,
    ) -> None:
        """
        Switch the global integrator to the user-specified one for n_steps steps.

        Use this by adding it to the list of reporters, typically near the end.
        It hooks itself to the on_reaction Simulation hook.

        :param minimizer: An OpenMM integrator, typically a LangevinMiddleIntegrator,
            with a high friction and smaller timestep, but in principle can be any.
        :param n_steps: How many extra steps to perform with the given integrator.
        """
        self.minimizer = minimizer
        self.n_steps = n_steps
        self.integrator_index = None

    def pre_simulation_start(self, simulation: Simulation) -> None:
        assert simulation.integrator is not None
        self.integrator_index = simulation.integrator.addIntegrator(self.minimizer)

    def on_reaction(
        self, simulation: Simulation, reactions: list[tuple[str, list[Fragment]]]
    ) -> None:
        prev_integrator = simulation.context.get_current_integrator()
        assert self.integrator_index is not None
        simulation.context.set_current_integrator(self.integrator_index)
        simulation.context.do_steps(self.n_steps)
        simulation.context.set_current_integrator(prev_integrator)
