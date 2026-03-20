class Reporter:
    def on_simulation_start(self, simulation):
        """
        Called once when Simulation is constructed.
        """
        pass

    def on_simulation_finish(self, simulation):
        """
        Called once when simulation's finish() is called.
        """
        pass

    def on_trajectory_frame(self, simulation):
        """
        Called every traj_frequency frames. There is a single per simulation traj_frequency because that's a simple
        way of getting multiple output types with nicely aligned time frames.
        """
        pass

    def interactive_line(self, simulation) -> str:
        """
        Should return its addition to the interactive status progress display.
        """
        pass

    def pre_modification(self, simulation):
        """
        Called before the modification algorithm, but only if there is any reactions happening.
        """

    def on_reaction(self, simulation, reactions):
        """
        Called after the modification algorithm runs. Reactions is the list of reactions that were applied.
        """
        pass

    def post_reaction(self, simulation):
        """
        Called after minimization (if done) and reinitialize.
        """
        pass