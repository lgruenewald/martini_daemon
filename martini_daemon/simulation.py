import openmm as mm
import os
import zlib
import sys
from datetime import datetime
from importlib.metadata import version
from typing import Any

from .__formats import *
from .__parser import GromacsTopFile, InvalidTopologyError
from .__core import System, Context, wrap_coupling


class Simulation:
    """
    Simulation class.

    * provide a friendly interface for reporters requesting output files. contains default values for Martini simulations.
    * hold simulation metadata, such as current step, simulation name.
    * own all open file handles, reporters and loggers.
    * be passed around to all reporters to provide the required metadata and file handle access for reporting.
    * is the required glue between all components, also only uses the public interface of different components.
    * provide access to the system, context, topstar instance to all reporters and the user
    """

    def __init__(
        self, top_path: str, geom_path: str, md_steps: int,
        reporters: list[Reporter],
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
    ):
        """

        """
        self.reporters = []
        self.step = 0
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


        # file handles setup
        # dict of suffix -> (handle, compression_obj | None)
        self.output_files = {}

        self.open(".log")
        self.info("Simulation __init__ called")
        self.info("Martini Daemon version", version("martini_daemon"))
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

        # Parsing
        self.info("Parsing start")
        box, start_pos, start_vel = read_geometry(geom_path)
        options["respos"] = read_geometry(restraint_coord_path)[1] if restraint_coord_path is not None else start_pos
        self.info(f"Read {len(start_pos)} atoms from {geom_path}. Box: {box.to_lattice()}. Velocities read? {start_vel is not None}.")
        self.system = System(options=options)
        try:
            top_parser = GromacsTopFile(self.system, top_path)
        except InvalidTopologyError:
            self.error("Fatal error during .top parsing.")
            # I want a silent exit, kinda hacky..
            raise SystemExit
        self.info("Parsing finished")

        self.info("setup integrator", "dt (ns):", self.dt_ns, "type:", type(self.md_integrator).__name__)
        # integrator -> always compound, index 0 always for md
        self.integrator = mm.CompoundIntegrator()
        self.integrator.addIntegrator(self.md_integrator)
        # couplings
        for c in coupling:
            self.system.add_force(wrap_coupling(c)(self.system))

        # build context
        self.info("Building context")
        self.context = Context(self.system, self.integrator, platform, context_parameters)

        # set pos, vel
        self.context.set_positions(start_pos, box)
        if start_vel is not None:
            self.context.set_velocities(start_vel)

        # TODO reporters


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
            f"Simulation {self.__sim_name}, step {self.step}, time {self.time_ns}.",
            self.system.get_atom_names(),
            self.system.get_res_ids(),
            self.system.get_res_names(),
            box, pos, vel
        )

        # TODO replay
        # TODO step
        # TODO simulate - make sure to call finish after
