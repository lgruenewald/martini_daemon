from .__rust import Fragment


class Reporter:
    def pre_simulation_start(self, simulation) -> None:
        """Called just before the context is initialized.

        Use on_simulation_start unless you really need to mutate simulation in a way that needs to happen
        before context initialization.
        """

    def on_simulation_start(self, simulation, continue_sim: bool) -> None:
        """Called once when Simulation is constructed. After the context is initialized.

        :param continue_sim: If True, should append instead of overwrite. Warning! May need to truncate files to
        simulation.current_step first! Do not write headers twice! Prefer to raise NotImplementedError if truncating
        is needed, but it is not implemented.
        """
        pass

    def on_simulation_finish(self, simulation) -> None:
        """Called once when simulation's finish() is called."""
        pass

    def on_trajectory_frame(self, simulation) -> None:
        """Called every traj_frequency frames. There is a single per simulation traj_frequency because that's a simple
        way of getting multiple output types with nicely aligned time frames.
        """
        pass

    def interactive_line(self, simulation) -> str | None:
        """Should return its addition to the interactive status progress display."""
        pass

    def pre_modification(self, simulation) -> None:
        """Called before the modification algorithm, but only if there may be any reactions happening."""

    def on_reaction(
        self, simulation, reactions: list[tuple[str, list[Fragment]]]
    ) -> None:
        """Called after the modification algorithm runs. Reactions is the list of reactions that were applied."""
        pass

    def post_reaction(self, simulation) -> None:
        """Called after all on_reaction reporters were resolved (some might apply minimization)."""
        pass
