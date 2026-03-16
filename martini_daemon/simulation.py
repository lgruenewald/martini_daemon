import openmm as mm
import os
import zlib
import sys
from datetime import datetime
from importlib.metadata import version
from .__formats import *


class Simulation:
    """
    Simulation class.

    * hold simulation metadata, such as current step, simulation name.
    * own all open file handles and loggers.
    * provide a friendly interface for reporters requesting output files.
    * be passed around to all reporters to provide the required metadata for reporting.
    * own all reporters, and provide an API for calling them all whenever a certain event occurs.
    * provide access to the system, context, topstar instance to all reporters and the user
    """

    def __init__(
        self, top_path: str, geom_path: str, md_steps: int,
        dm_frequency: int = 0, traj_frequency: int = 0,
        sim_name: str = "out", continue_sim=False,
        reporters = None,
        coupling = None,
        integrator: mm.Integrator | None = None,
        epsilon_r: float = 15., cutoff_nm: float = 1.1,
        include_dir: str = None,
        defines: dict[str, str] = None,
        platform: str | None | mm.Platform = None,
        context_parameters: None | dict[str, str] = None,
        restraint_coord_path=None,
    ):
        # TODO continue_sim implementation
        self.reporters = []
        self.step = 0
        self.total_steps = md_steps
        self.sim_name = sim_name
        self.time_ns = 0.
        self.reporters = reporters or []
        self.integrator = integrator or mm.LangevinMiddleIntegrator(
            300 * mm.unit.kelvin, 1. / mm.unit.picosecond, 0.02 * mm.unit.picosecond
        )
        self.coupling = coupling
        if self.coupling is None:
            self.coupling = [
                mm.MonteCarloBarostat(
                    1.0 * mm.unit.bar,
                    300 * mm.unit.kelvin
                ),
                mm.CMMotionRemover()
            ]
        self.dt_ns: float = (
            self.integrator.getStepSize().value_in_unit(mm.unit.nanosecond)
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

        _, res_pos, _ = read_geometry(restraint_coord_path or geom_path)

        # file handles setup
        # dict of suffix -> (handle, compression_obj | None)
        self.output_files = {}

        self.open(".log")
        self.info("Simulation __init__ called")
        self.info("Martini Daemon version", version("martini_daemon"))
        self.info(
            "Parameters:", top_path, geom_path, "dt (ns):", self.dt_ns,
            "steps:", self.total_steps, "dm_freq:", self.dm_frequency,
            "traj_freq:", self.traj_frequency, "sim_name:", self.sim_name,
            "platform:", platform, "context_parameters:", context_parameters,
            "defines:", defines, "include_dir:", include_dir,
        )

        # Parsing
        self.info("Parsing start")


        self.info("Parsing finished")


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
        path = self.sim_name + suffix
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
        Flushes output files and closes output file handles.
        """
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

    # Friendly interface for requesting output files -- trajectory, geometry

    # Friendly interface for setting up and running simulations

