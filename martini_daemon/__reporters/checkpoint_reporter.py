import os

from ..__formats import write_checkpoint
from ..__reporter import Reporter
from ..__simulation import Simulation


class CheckpointReporter(Reporter):
    def __init__(self, n_checkpoints: int = 3) -> None:
        """Create a Checkpoint Reporter.

        The Checkpoint Reporter will write checkpoint files, allowing for continuable simulations in case of
        interruption. Load checkpoints by instantiating the class CheckpointLoader(). Using the same arguments as
        you would a simulation, but with the path to the .chk file as the geom_path argument.

        Warning! Read the following notes before using.

        * The best solution to checkpoints, with an acceptable performance, granularity over frequencies
            and a proper format for storing topologies, including user-defined forces is still under consideration.
        * Currently, the checkpoints are written every trajectory frame and only store atom positions, velocities,
            simulation frame and time.
        * Reporter output files will be truncated to the appropriate frame at which the checkpoint was taken.
            A backup copy is made before truncation. After this, reporters will append their output to the truncated
            output files.
        * The simulation topology is recovered based on re-parsing the .top file, and replaying all reactions based on
            the output of ReactionReporter. If you intend to continue simulations with reactions, you MUST also have
            ReactionReporter included in the list of reporters.
        * Only the checkpoint files are atomic!
            First the old checkpoints are renamed to increment the number in their name.
            The newest .chk file is only renamed to have the extension .chk if it is complete, otherwise it is .chk_tmp.
            Other reporters output may not be atomic. To increase the chances of having good trajectory outputs
            in case of sudden power failure, the CheckpointReporter should be the last reporter provided to
            Simulation, so the new checkpoint is only written if all other reporters finish writing a trajectory frame.
            These checkpoints are on a best-effort basis, they should be only used as a last resort.
            Ideally, the user should also first verify the state of the trajectories before continuing a simulation,
            and make manual backups.
        * Simulation settings provided to the checkpoint loading function must be identical to the arguments
            to the original Simulation() constructor.
        * Reloading is not deterministic. Random state for e.g. coupling schemes involving randomness will be lost.
        * The format of checkpoint files can change between minor versions without warning. Only load checkpoints
            with the same version as they were written with. Thanks to a magic number at the beginning of the file,
            loading a checkpoint file with the wrong version should raise an appropriate error to the user.

        :param n_checkpoints: Number of checkpoint files to keep at a time. '.chk' is the newest,
            '.chk2' the one before, with increasing numbers representing older checkpoint files.
        """
        self.n_checkpoints = n_checkpoints
        self.paths: list[str] = []
        self.tmp_path: str = ""

    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        self.paths = [
            simulation.request_path(f".chk{i + 1 if i > 0 else ''}", copy=True)
            for i in range(self.n_checkpoints + 1)
        ]
        self.tmp_path: str = simulation.request_path(".chk_tmp")

    def on_trajectory(self, sim: Simulation) -> None:
        # increment the number on checkpoint files
        for i in range(self.n_checkpoints, 0, -1):
            # if n=3, will count down as 3, 2, 1
            os.rename(self.paths[i - 1], self.paths[i])

        # make the new checkpoint
        pos, box = sim.context.get_positions()
        vel = sim.context.get_velocities()
        write_checkpoint(
            self.tmp_path,
            sim.current_step,
            sim.trajectory_frame,
            sim.reactions_so_far,
            sim.time_ps,
            len(pos),
            box,
            pos,
            vel,
        )

        # commit to the new checkpoint atomically by renaming it
        os.rename(self.tmp_path, self.paths[0])

        # delete the oldest checkpoint only as a final step
        os.remove(self.paths[-1])


class CheckpointLoader(Simulation):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


# TODO - check if .chk is provided, warn if there are reactions in the .chk but not continue_sim
