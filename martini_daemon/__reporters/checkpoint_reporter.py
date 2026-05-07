import os
import sys

from ..__formats import read_checkpoint, write_checkpoint
from ..__reporter import Reporter
from ..__simulation import Simulation
from .reaction_reporter import ReactionReporter


class CheckpointReporter(Reporter):
    def on_simulation_finish(self, simulation) -> None:
        pass

    def __init__(self, interval: int = 0, n_checkpoints: int = 3) -> None:
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

        :param interval: Lower bound on the number of MD steps to wait between checkpoints. Will only make checkpoints
            during trajectory frames if the number of frames since last checkpoint is larger than this.
        :param n_checkpoints: Number of checkpoint files to keep at a time. '.chk' is the newest,
            '.chk2' the one before, with increasing numbers representing older checkpoint files.

        """
        self.n_checkpoints = n_checkpoints
        self.paths: list[str] = []
        self.tmp_path: str = ""
        self.interval = interval
        self.__last_chk = 0

    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        self.paths = [
            simulation.request_path(f".chk{i + 1 if i > 0 else ''}", continue_sim=continue_sim)
            for i in range(self.n_checkpoints + 1)
        ]
        self.tmp_path: str = simulation.request_path(".chk_tmp")
        # no chk at the start
        self.__last_chk = simulation.current_step

    def on_trajectory_frame(self, sim: Simulation) -> None:
        if sim.current_step < self.__last_chk + self.interval:
            # interval condition not fulfilled yet
            return
        self.__last_chk = sim.current_step
        # increment the number on checkpoint files
        for i in range(self.n_checkpoints, 0, -1):
            # if n=3, will count down as 3, 2, 1
            if os.path.exists(self.paths[i - 1]):
                os.rename(self.paths[i - 1], self.paths[i])

        # make the new checkpoint
        pos, box = sim.context.get_positions()
        vel = sim.context.get_velocities()
        write_checkpoint(
            self.tmp_path,
            # hacky, but simulation does not expose sim name otherwise so reporters are forced to use request_path
            self.tmp_path.removesuffix(".chk_tmp"),
            sim.current_step,
            # the index of current trajectory frame being written
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
        if os.path.exists(self.paths[-1]):
            os.remove(self.paths[-1])


class CheckpointLoader(Simulation):
    def __init__(self, checkpoint: str, *args, **kwargs) -> None:
        """Load a checkpoint from chk_path.

        Pass additional arguments as you would to Simulation().
        """
        chk = read_checkpoint(checkpoint)
        sim_name = chk.sim_name
        reactions = None
        if chk.reactions_so_far > 0:
            rx_path = sim_name + ".reactions"
            assert os.path.exists(rx_path), (
                f"{rx_path} does not exist, but reactions need to be replayed."
            )
            reactions = ReactionReporter.read_reactions(rx_path)
        assert "checkpoint" not in kwargs, "Use chk_path, not checkpoint."
        super().__init__(*args, checkpoint=chk, **kwargs)
        # REPLAY
        print()
        if reactions is not None:
            for step, rx_name, frags in reactions:
                if self.current_step < step:
                    break
                frag_ids = [frag_id for (name, frag_id, atoms) in frags]
                sys.stdout.write(f"\033[2K\rReplaying reactions: {step}/{self.current_step}: {rx_name} {frag_ids}")
                # verification
                frag_objs = [
                    self.top.frag_list.get_fragment(frag_id) for frag_id in frag_ids
                ]
                for (name, frag_id, atoms), frag_obj in zip(frags, frag_objs):
                    assert name == frag_obj.name
                    assert frag_id == frag_obj.frag_id
                    assert all(a == b for a, b in zip(atoms, frag_obj.atoms))
                self.top.modification([(rx_name, frag_ids)])
            sys.stdout.write(f"\033[2K\rReplaying reactions: done replaying {len(reactions)} reactions.")
        print()
