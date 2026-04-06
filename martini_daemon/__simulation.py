import openmm as mm
import os
import zlib
import math
import sys
from datetime import datetime
from time import time
from importlib.metadata import version
from typing import Any, Type, Callable

from .__formats import *
from .__parser import GromacsTopFile, InvalidTopologyError
from .__core import System, Context, wrap_coupling
from .__forces import NonBonded
from .__rust import build_version, Fragment
from .__reporter import Reporter
from .__topstar import TopStar


class Simulation:
    """
    Simulation class.

    * provide a friendly interface for reporters requesting output files. contains default values for Martini simulations.
    * hold simulation metadata, such as current step, simulation name.
    * own all open file handles, reporters and loggers.
    * be passed around to all reporters to provide the required metadata and file handle access for reporting.
    * is the required glue between all components, also only uses the public interface of different components.
    * provide access to the system, context, topstar instance to all reporters and the user

    geom_path: path to geometry. Note: if None is passed, no context will be initialized. Simulation can then be
    used as a .top parser, and the resulting topology can then be read out. Of course things needing Context
    steps will not be available.
    """

    def __init__(
        self, top_path: str, geom_path: str | None, md_steps: int,
        reporters: list[Reporter] | None = None,
        dm_frequency: int = 0, traj_frequency: int = 0,
        sim_name: str = "out", continue_sim=False,
        coupling = None,
        integrator: mm.Integrator | None = None,
        options: dict[str, Any] | None = None,
        include_dir: str = None,
        defines: dict[str, str] = None,
        platform: str | None | mm.Platform = None,
        context_parameters: None | dict[str, str] = None,
        restraint_coord_path=None,
        nonbonded: Callable[[], NonBonded] | Type[NonBonded] | None = None,
    ):
        """

        """
        self.reporters = reporters or []
        self.current_step = 0
        self.total_steps = md_steps
        self.__sim_name = sim_name
        self.time_ns = 0.
        self.md_integrator = integrator or mm.LangevinMiddleIntegrator(
            300 * mm.unit.kelvin, 1. / mm.unit.picosecond, 0.02 * mm.unit.picosecond
        )
        if coupling is None:
            coupling = [
                mm.MonteCarloBarostat(
                    1.0 * mm.unit.bar,
                    300 * mm.unit.kelvin
                ),
                mm.CMMotionRemover()
            ]
        self.dt_ns: float = (
            self.md_integrator.getStepSize().value_in_unit(mm.unit.nanosecond)
        )
        self.dm_frequency = dm_frequency
        self.traj_frequency = traj_frequency
        if type(platform) is str:
            platform = mm.Platform.getPlatformByName(platform)
        include_dir = include_dir or (
            "GMXDATA" in os.environ and
            os.path.join(os.environ["GMXDATA"], "top")
        ) or (
            "GMXBIN" in os.environ and
            os.path.join(os.environ["GMXBIN"], "..", "share", "gromacs", "top")
        ) or "/usr/local/gromacs/share/gromacs/top"
        if nonbonded is None:
            nonbonded = NonBonded

        # file handles setup
        # dict of suffix -> (handle, compression_obj | None)
        self.output_files = {}

        self.open(".log")
        self.info(f"Martini Daemon {version('martini_daemon')} log file")
        self.info("Build version:", build_version())
        self.info(
            "Parameters:", top_path, geom_path,
            "steps:", self.total_steps, "dm_freq:", self.dm_frequency,
            "traj_freq:", self.traj_frequency, "sim_name:", self.__sim_name,
            "platform:", platform, "context_parameters:", context_parameters,
            "defines:", defines, "include_dir:", include_dir,
        )
        if options is None:
            options = {}
        if options.get("epsilon_r") is None:
            options["epsilon_r"] = 15.
        if options.get("cutoff") is None:
            options["cutoff"] = 1.1
        if options.get("respos") is None:
            if restraint_coord_path is not None:
                options["respos"] = read_geometry(restraint_coord_path)[1]
            elif geom_path is not None:
                options["respos"] = read_geometry(geom_path)[1]

        # Parsing
        self.info("Parsing start")
        self.system = System(options=options)
        try:
            GromacsTopFile(self.system, top_path)
        except InvalidTopologyError:
            self.error("Fatal error during .top parsing.")
            # I want a silent exit, kinda hacky..
            raise SystemExit
        nb = nonbonded(self.system)
        excl = nb.get_exclusion_helper()
        self.system.add_force(nb)
        self.system.add_force(excl)
        self.system.build_initial_molecules()
        self.info("Parsing finished")

        self.info("TopStar build start")
        self.top = TopStar(self.system)
        self.info("TopStar build finished")

        self.info("setup integrator", "dt (ns):", self.dt_ns, "type:", type(self.md_integrator).__name__)
        # integrator -> always compound, index 0 always for md
        self.integrator = mm.CompoundIntegrator()
        self.integrator.addIntegrator(self.md_integrator)
        # couplings
        for c in coupling:
            self.system.add_force(wrap_coupling(c)(self.system))

        # build context
        if geom_path is not None:
            box, start_pos, start_vel = read_geometry(geom_path)
            self.info(f"Read {len(start_pos)} atoms from {geom_path}. Box: {box.to_lattice()}. Velocities read? {start_vel is not None}.")

            self.info("Building context")
            self.context = Context(self.system, self.integrator, box, platform, context_parameters)

            # set pos, vel
            self.context.set_positions(start_pos, box)
            if start_vel is not None:
                self.context.set_velocities(start_vel)

        for r in self.reporters:
            r.on_simulation_start(self)

        # for the estimated time left display
        self.last_step_time = 0.
        self.first_step_time = 0.
        self.trajectory_frame = 0

    # File handles and loggers
    @staticmethod
    def __backup_try(path):
        parent, filename = os.path.split(path)
        if os.path.isfile(path):
            bkup_num = 0
            bkup_path = path
            while os.path.isfile(bkup_path):
                bkup_num += 1
                bkup_path = os.path.join(parent, f"#{filename}.{bkup_num}#")
            os.rename(path, bkup_path)
            print(f"Backed up {path} to {bkup_path}")

    def request_path(self, suffix):
        """
        Convert suffix to path.
        """
        path = self.__sim_name + suffix
        self.__backup_try(path)
        return path

    def open(self, suffix, compress=False):
        """
        Opens a new file handle for writing.
        """
        path = self.request_path(suffix)
        self.output_files[suffix] = (open(path, "wb"), zlib.compressobj(6) if compress else None)

    def write(self, suffix, bytes_or_text):
        """
        Writes to open handle.
        """
        if type(bytes_or_text) is str:
            bytes_or_text = bytes_or_text.encode("utf-8")
        if self.output_files[suffix][1] is not None:
            bytes_or_text = self.output_files[suffix][1].compress(bytes_or_text)
        self.output_files[suffix][0].write(bytes_or_text)

    def print(self, suffix, *args, sep=" ", end="\n"):
        """
        Writes all args to output file, separated by separator (space). Writes a newline after.
        """
        for i, arg in enumerate(args):
            if i > 0:
                self.write(suffix, sep)
            self.write(suffix, f"{arg}")
        self.write(suffix, end)

    def flush(self):
        """
        Flushes all output files. Also called at trajectory frames automatically.
        """
        for handle, comp in self.output_files.values():
            if comp is not None:
                handle.write(comp.flush())
            handle.flush()

    def finish(self):
        """
        1. Runs finish on all reporters.
        2. Flushes output files and closes output file handles.
        """
        for r in self.reporters:
            r.on_simulation_finish(self)
        self.flush()
        for handle, _ in self.output_files.values():
            handle.close()
        self.output_files = {}

    def info(self, *args):
        self.print(
            ".log",
            datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S,%f")[:-3],
            *args
        )

    def warn(self, message):
        self.info(
            "[WARNING] " + message
        )
        print("[WARNING]", message, file=sys.stderr)


    def error(self, message):
        self.info(
            "[ERROR] " + message
        )
        print("[ERROR]", message, file=sys.stderr)


    # Friendly interface for setting up and running simulations
    @staticmethod
    def set_process_title(newname=b"daemon"):
        """
        Set process title to something else than "python".
        Nothing critical, just a nice thing to have for top/htop/btop/mu ;)
        Only works on (some versions of) linux.
        May fail silently with no exceptions thrown.

        based on: https://stackoverflow.com/questions/564695/is-there-a-way-to-change-effective-process-name-in-python

        :param newname: new process name, as a byte string.
        """

        try:
            from ctypes import cdll, byref, create_string_buffer
            libc = cdll.LoadLibrary('libc.so.6')
            buff = create_string_buffer(len(newname)+1)
            buff.value = newname
            libc.prctl(15, byref(buff), 0, 0, 0)
        finally:
            pass

    def save_geometry(self, path) -> None:
        """
        Save the current simulation geometry to path.

        This may include atom names, residue id, residue names, timestep, current time,
        simulation name, positions, velocities and the pbc box, depending on the file format used.

        :param path: path to save geometry to.
        """
        pos, box = self.context.get_positions()
        vel = self.context.get_velocities()
        write_geometry(
            path,
            f"Simulation {self.__sim_name}, step {self.current_step}, time {self.time_ns}.",
            self.system.get_atom_names(),
            self.system.get_res_names(),
            self.system.get_res_ids(),
            box, pos, vel
        )

        # TODO replay
        # TODO checkpoints

    def simulate(self):
        remaining = self.total_steps - self.current_step
        sim_ns = self.total_steps * self.dt_ns
        print(f"Simulation of {remaining} steps ({sim_ns} ns)")
        gcd = math.gcd(self.traj_frequency, self.dm_frequency, remaining)
        print(f"D/M freq {self.dm_frequency} Traj freq {self.traj_frequency} gcd {gcd}")
        while self.current_step < self.total_steps:
            self.step(
                gcd, traj=self.current_step % self.traj_frequency == 0 if self.traj_frequency > 0 else False,
                dm=self.current_step % self.dm_frequency == 0 if self.dm_frequency > 0 else False
            )
        self.do_traj_frame()
        print()
        self.finish()

    @staticmethod
    def __format_time(total):
        c = int(total)
        seconds = c % 60
        c //= 60
        minutes = c % 60
        c //= 60
        hours = c % 24
        c //= 24
        days = c

        res = f"{hours:02}:{minutes:02}:{seconds:02}"
        if days > 0:
            res = f"{days}d " + res
        return res

    def do_traj_frame(self):
        for r in self.reporters:
            self.info(f"Trajectory {r.__class__.__name__} start")
            r.on_trajectory_frame(self)
            self.info(f"Trajectory {r.__class__.__name__} finished")
        self.trajectory_frame += 1

    def step(self, n_steps: int, traj=False, dm=False, silent=False):
        """
        Do the following:
        - steps n_steps
        - D/M algorithm if dm is true
        - locally minimize energy and reinitialize system if reactions happened
        - trigger a trajectory frame on reporters if traj is True
        - display info to logs and screen, % info given by self.current_step and self.md_steps
        - update self.current_step
        """
        start_time = time()
        if traj:
            self.do_traj_frame()

        self.info(f"doing md steps to go from {self.current_step} to")
        self.current_step += n_steps
        percent = self.current_step / self.total_steps * 100. if self.total_steps > 0 else 100.
        self.info(f"step {self.current_step}")
        if n_steps > 0:
            self.info(f"md_steps {n_steps}")
            self.info("MD start")
            try:
                self.context.do_steps(n_steps)
            except mm.OpenMMException as e:
                self.error(f"!!! OpenMM Exception !!!\n{e}")
            self.info("MD finished")
        if self.total_steps > 0 and n_steps > 0 and not silent:
            self.time_ns = self.dt_ns * self.current_step
            time_left = self.__format_time(self.last_step_time * (self.total_steps - self.current_step))
            reporter_data = " ".join(filter(None, [r.interactive_line(self) for r in self.reporters]))
            sys.stdout.write(
                f"\033[2K\rstep {self.current_step}"
                f"({self.time_ns:.2f} ns, "
                f"{percent:.1f}%) "
                f"{time_left} {reporter_data}"
            )
        if dm:
            self.info("Detection start")
            pos, box = self.context.get_positions()
            reactions: list[tuple[str, list[int]]] = self.top.detection(box, pos)
            self.info(f"After detection there were {len(reactions)} reactions")
            self.info("Detection finished")
            if len(reactions) > 0:
                for r in self.reporters:
                    r.pre_modification(self)
                # we now need copies of fragments, since they possibly got consumed in the reaction
                reactions: list[tuple[str, list[Fragment]]] = self.top.modification(reactions)
                for r in self.reporters:
                    r.on_reaction(self, reactions)
                # TODO minimization
        end_time = time()
        if n_steps > 0:
            step_time = (end_time - start_time) / n_steps
            if self.last_step_time > 0.:
                self.last_step_time = step_time * 0.01 + self.last_step_time * 0.99
            elif not traj or not dm:
                # step 0 tends to have both xtc and dm as True, and is
                # usually unrepresentatively slow
                scale = self.dm_frequency / self.traj_frequency
                if scale > 1.:
                    scale = 1. / scale
                if self.first_step_time == 0.:
                    # continuations might not start with an expensive step
                    scale = 0.
                self.last_step_time = step_time * (1 - scale) + scale * self.first_step_time
            else:
                # step 0 probably
                self.first_step_time = step_time
