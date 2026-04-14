import math
import os
import shutil
import sys
import traceback
from collections.abc import Callable
from datetime import datetime
from importlib.metadata import version
from time import time
from typing import Any, BinaryIO

import openmm as mm
import openmm.app as mmapp
from openmm.unit import md_unit_system

from .__core import Context, System, wrap_coupling
from .__forces import NonBonded
from .__formats import read_geometry, write_geometry
from .__parser import GromacsTopFile, InvalidTopologyError
from .__reporter import Reporter
from .__rust import Fragment, build_version
from .__topstar import TopStar


class Simulation:
    def __init__(
        self,
        top_path: str,
        geom_path: str | None,
        md_steps: int,
        reporters: list[Reporter] | None = None,
        dm_frequency: int = 0,
        traj_frequency: int = 0,
        sim_name: str = "out",
        coupling=None,
        integrator: mm.Integrator | None = None,
        options: dict[str, Any] | None = None,
        include_dirs: list[str] | None = None,
        defines: dict[str, str] | None = None,
        platform: str | None | mm.Platform = None,
        context_parameters: None | dict[str, str] = None,
        nonbonded: Callable[[System], NonBonded] | type[NonBonded] | None = None,
        continue_sim: bool = False,
    ):
        """Simulation class.

        * Provides a friendly interface for reporters requesting output files. Contains default values for Martini simulations.
        * Holds simulation metadata, such as current step, simulation name.
        * Owns all open file handles, reporters and loggers.
        * Is passed around to all reporters to provide the required metadata and file handle access for reporting.
        * Is the required glue between all components, also only uses the public interface of different components.
        * Provides access to the system, context, topstar instance to all reporters and the user

        :param top_path: Path to the Martini .top file.
        :param geom_path: Path to the geometry file (.gro, .xyz).
        :param md_steps: Number of MD steps.
        :param reporters: List of reporters to use during simulation.
        :param dm_frequency: Frequency of the Detection/Modification algorithm.
        :param traj_frequency: Frequency of Trajectory frames.
        :param sim_name: Short name of the simulation. All output files will be prefixed by this name.
        :param coupling: List of OpenMM coupling forces to use. If None, pressure coupling at 1 bar and 300 kelvin,
            and center of mass motion removal will be employed.
        :param integrator: Base integrator to use during the simulation. Note: a compound integrator will be set up
            based on it. Depending on the reporters, local minimization or other integrators can be configured alongside.
            If None, LangevinMiddleIntegrator will be used, at 300 kelvin, 1 ps-1 collision frequency and 0.02 ps dt.
        :param options: Additional data to pass to the system. Example keys available are "epsilon_r" (default 15),
            "cutoff" (default 1.1, in nanometers) to control the nonbonded force, as well as "respos",
            which can be set as a f64 (n_atoms, 3) shaped numpy array for position restraint reference coordinates
            (default same as geom_path coordinates).
        :param include_dirs: Additional include directories for #include directives in .top files. By default it tries
            to detect the gromacs installation and add an entry to the "top" subfolder inside it.
        :param defines: Additional defines to pass to the .top parser.
        :param platform: Which OpenMM platform to use.
        :param context_parameters: Additional options to pass to the platform.
        :param nonbonded: Nonbonded force to use, passed as a type or a function that returns the martini daemon Force
            when called with system as its argument. By default, the Martini compatible shifted Lennard-Jones
            and reaction-field electrostatics are used.
        :param continue_sim: Whether this is a continuation of a previous simulation. Do not use manually! Use
            ReactionReporter's replay method to continue simulations.

        Note: You may want to take a look at the following attributes, which also contain methods for common simulation
        tasks:

        :attribute context: See :doc:`/autoapi/martini_daemon/Context`. Note: if geom_path is None, no context will be initialized.
            Simulation can then be solely used as a .top parser, the resulting topology and OpenMM system can still be read out.
        :attribute system: See :doc:`/autoapi/martini_daemon/System`.
        """
        self.__reporters = reporters or []
        self.current_step: int = 0
        self.total_steps: int = md_steps
        self.__sim_name = sim_name
        self.time_ps: float = 0.0
        md_integrator = integrator or mm.LangevinMiddleIntegrator(
            298.0,  # kelvin
            1.0,  # ps^-1
            0.02,  # ps
        )
        if coupling is None:
            coupling = [
                mm.MonteCarloBarostat(
                    1.0,  # bar
                    298.0,  # kelvin
                ),
                mm.CMMotionRemover(),
            ]
        self.dt_ps: float = md_integrator.getStepSize().value_in_unit_system(
            md_unit_system
        )  # ps
        self.dm_frequency: int = dm_frequency
        self.traj_frequency: int = traj_frequency
        if type(platform) is str:
            platform = mm.Platform.getPlatformByName(platform)
        if include_dirs is None:
            include_dirs = (
                (
                    "GMXDATA" in os.environ
                    and [os.path.join(os.environ["GMXDATA"], "top")]
                )
                or (
                    "GMXBIN" in os.environ
                    and [
                        os.path.join(
                            os.environ["GMXBIN"], "..", "share", "gromacs", "top"
                        )
                    ]
                )
                or ["/usr/local/gromacs/share/gromacs/top"]
            )
        if nonbonded is None:
            nonbonded = NonBonded

        # file handles setup
        # dict of suffix -> (handle, compression_obj | None)
        self.__output_files: dict[str, BinaryIO] = {}

        self.open(".log")
        self.info(f"Martini Daemon {version('martini_daemon')} log file")
        self.info("Build version:", build_version())
        self.info(
            "Parameters:",
            top_path,
            geom_path,
            "steps:",
            self.total_steps,
            "dm_freq:",
            self.dm_frequency,
            "traj_freq:",
            self.traj_frequency,
            "sim_name:",
            self.__sim_name,
            "platform:",
            platform,
            "context_parameters:",
            context_parameters,
            "defines:",
            defines,
            "include_dirs:",
            include_dirs,
        )
        if options is None:
            options = {}
        if options.get("epsilon_r") is None:
            options["epsilon_r"] = 15.0
        if options.get("cutoff") is None:
            options["cutoff"] = 1.1
        if options.get("respos") is None and geom_path is not None:
            options["respos"] = read_geometry(geom_path)[1]

        # Parsing
        self.info("Parsing start")
        self.system: System = System(options=options)
        try:
            GromacsTopFile(
                self.system, top_path, include_dirs=include_dirs, defines=defines
            )
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
        self.top: TopStar = TopStar(self.system)
        self.info("TopStar build finished")

        self.info(
            "setup integrator",
            "dt (ps):",
            self.dt_ps,
            "type:",
            type(md_integrator).__name__,
        )
        self.integrator: mm.CompoundIntegrator | None = (
            mm.CompoundIntegrator()
        )  #: Compound Integrator with integrator index 0 as the user specified integrator. Only exposed before the context is built.
        self.integrator.addIntegrator(md_integrator)
        # couplings
        for c in coupling:
            self.system.add_force(wrap_coupling(c)(self.system))

        # Reporters!
        # reporters can add integrators only here
        for r in self.__reporters:
            r.pre_simulation_start(self)

        self.__context: Context | None = None
        # build context
        if geom_path is not None:
            box, start_pos, start_vel = read_geometry(geom_path)
            self.info(
                f"Read {len(start_pos)} atoms from {geom_path}. Box: {box.to_lattice()}. Velocities read? {start_vel is not None}."
            )

            self.info("Building context")
            assert platform is None or isinstance(platform, mm.Platform)
            self.__context = Context(
                self.system, self.integrator, box, platform, context_parameters
            )

            # set pos, vel
            self.__context.set_positions(start_pos, box)
            if start_vel is not None:
                self.__context.set_velocities(start_vel)

        # it's owned by context now
        self.integrator = None

        for r in self.__reporters:
            r.on_simulation_start(self, continue_sim)

        # for the estimated time left display
        self.__last_step_time = 0.0
        self.__first_step_time = 0.0
        self.trajectory_frame: int = 0

        # to avoid double finish
        self.__finished = False

    # File handles and loggers
    @staticmethod
    def __backup_try(path: str, copy: bool = False) -> None:
        parent, filename = os.path.split(path)
        if os.path.isfile(path):
            bkup_num = 0
            bkup_path = path
            while os.path.isfile(bkup_path):
                bkup_num += 1
                bkup_path = os.path.join(parent, f"#{filename}.{bkup_num}#")
            if copy:
                os.rename(path, bkup_path)
            else:
                shutil.copy(path, bkup_path)
            print(f"Backed up {path} to {bkup_path}")

    def request_path(self, suffix: str, copy: bool = False) -> str:
        """Request a writable path for an output file. Back up the file if it already exists.

        :param suffix: suffix to use. Usually a file extension, e.g. ".xtc".
        :param copy: If True, it will make a copy of the original contents at the original path.
            Generally, only pass False, if you intend to append to the contents.
        """
        path = self.__sim_name + suffix
        self.__backup_try(path, copy=copy)
        return path

    def open(self, suffix: str, append: bool = False) -> None:
        """Opens a new file handle for writing, and ties the file handle's lifetime to the Simulation object.

        Note, the true path of the file will be a combination of simulation name and suffix.

        Also note that if the file already exists, it will be first backed up, and the new path will be printed to
        stdout.

        Writing to these files should happen using Simulation.write() and Simulation.print().

        Simulation.finish() will automatically close it when the simulation ends. Simulation.close() should
        be called if an output file handle is no longer necessary before the simulation ends.

        If the handle is already open, it will only seek to the end, so that appending may resume.

        :param suffix: suffix to append to sim_name to get the path.
        :param append: if True, open file for appending, don't erase the file first.
        """
        if suffix in self.__output_files:
            self.__output_files[suffix].seek(0, os.SEEK_END)
        else:
            # if appending, make a copy instead
            path = self.request_path(suffix, copy=append)
            mode = "ab" if append else "wb"
            self.__output_files[suffix] = open(path, mode)  # noqa: SIM115

    def get_handle(self, suffix: str) -> BinaryIO:
        """Get file handle to do any operations on.

        Note: do not close the handle manually.
        Prefer Simulation.flush(), Simulation.close(), Simulation.write(), Simulation.print() for common operations.

        Note, the mode is always a binary mode. encode/decode manually from/to strings.
        """
        return self.__output_files[suffix]

    def write(self, suffix: str, bytes_or_text) -> None:
        """Writes to open handle."""
        if type(bytes_or_text) is str:
            bytes_or_text = bytes_or_text.encode("utf-8")
        assert self.__output_files.get(suffix) is not None
        self.__output_files[suffix].write(bytes_or_text)

    def print(self, suffix: str, *args, sep=" ", end="\n") -> None:
        """Writes all args to output file, separated by separator (space). Writes a newline after."""
        for i, arg in enumerate(args):
            if i > 0:
                self.write(suffix, sep)
            self.write(suffix, f"{arg}")
        self.write(suffix, end)
        self.flush(suffix)

    def flush(self, suffix: str) -> None:
        """Flushes a single output file handle."""
        assert self.__output_files.get(suffix) is not None
        self.__output_files[suffix].flush()

    def flush_all(self) -> None:
        """Flushes all output files. Also called at trajectory frames automatically."""
        for suffix in self.__output_files:
            self.flush(suffix)

    def close(self, suffix) -> None:
        """Closes a single open file."""
        self.flush(suffix)
        self.__output_files[suffix].close()
        del self.__output_files[suffix]

    def finish(self) -> None:
        """1. Runs finish on all reporters.
        2. Flushes output files and closes output file handles.

        Called by simulate(), or should be called by the user manually otherwise.
        """
        if self.__finished:
            return
        self.__finished = True
        for r in self.__reporters:
            r.on_simulation_finish(self)
        self.flush_all()
        for suffix in self.__output_files:
            self.close(suffix)

    def info(self, *args) -> None:
        """Writes a message to log."""
        self.print(
            ".log",
            datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S,%f")[:-3],
            *args,
        )

    def warn(self, message) -> None:
        """Writes a warning to the log file and to stderr."""
        self.info("[WARNING] " + message)
        print("[WARNING]", message, file=sys.stderr)

    def error(self, message) -> None:
        """Writes an error to the log file and stderr."""
        self.info("[ERROR] " + message)
        print("[ERROR]", message, file=sys.stderr)

    # Friendly interface for setting up and running simulations
    @staticmethod
    def set_process_title(newname=b"daemon") -> None:
        """Set process title to something else than "python". Nothing critical, purely aesthetic.
        Only works on (some versions of) linux. May fail silently with no exceptions thrown.

        Based on: https://stackoverflow.com/questions/564695/is-there-a-way-to-change-effective-process-name-in-python

        :param newname: new process name, as a byte string.
        """
        try:
            from ctypes import byref, cdll, create_string_buffer

            libc = cdll.LoadLibrary("libc.so.6")
            buff = create_string_buffer(len(newname) + 1)
            buff.value = newname
            libc.prctl(15, byref(buff), 0, 0, 0)
        finally:
            pass

    @property
    def context(self) -> Context:
        if self.__context is None:
            raise ValueError(
                "This Simulation has no context. Was Simulation() constructed with no geom_path?"
            )
        return self.__context

    def save_geometry(self, path) -> None:
        """Save the current simulation geometry to path.

        This may include atom names, residue id, residue names, timestep, current time,
        simulation name, positions, velocities and the pbc box, depending on the file format used.

        Supported file types include .gro and .xtc, this is determined based on the extension of the specified path.

        :param path: path to save geometry to.
        """
        pos, box = self.context.get_positions()
        vel = self.context.get_velocities()
        write_geometry(
            path,
            f"Simulation {self.__sim_name}, step {self.current_step}, time {self.time_ps} ps.",
            self.system.get_atom_names(),
            self.system.get_res_names(),
            self.system.get_res_ids(),
            box,
            pos,
            vel,
        )

        # TODO replay
        # TODO checkpoints

    def simulate(self, finish=True):
        """Performs the remaining steps (self.total_steps - self.current_step), and then calls self.finish() if finish is True.

        Will perform Detection/Modification and Trajectory writing according to their frequencies specified in
        the constructor for Simulation. Catches uncaught exceptions and logs them,
        and finishes the simulation prematurely if one occurs.

        :param finish: whether to close all output files. If False, must call finish() manually.
        """
        remaining = self.total_steps - self.current_step
        sim_ps = self.total_steps * self.dt_ps
        print(f"Simulation of {remaining} steps ({self.__format_sim_time(sim_ps)})")
        gcd = math.gcd(self.traj_frequency, self.dm_frequency, remaining)
        print(f"D/M freq {self.dm_frequency} Traj freq {self.traj_frequency} gcd {gcd}")
        try:
            while self.current_step < self.total_steps:
                self.step(
                    gcd,
                    traj=(
                        self.current_step % self.traj_frequency == 0
                        if self.traj_frequency > 0
                        else False
                    ),
                    dm=(
                        self.current_step % self.dm_frequency == 0
                        if self.dm_frequency > 0
                        else False
                    ),
                )
            self.__do_traj_frame()
            print()
        except Exception as e:
            traceback.print_exc()
            self.error(f"!!! Unexpected Exception!!!\n{e}")
        finally:
            if finish:
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

    @staticmethod
    def __format_sim_time(ps):
        if ps < 1000.0:
            return f"{ps:.2f} ps"
        if ps < 1000000.0:
            return f"{ps / 1000.0:.2f} ns"
        return f"{ps / 1000000:.2f} μs"

    def __do_traj_frame(self):
        for r in self.__reporters:
            self.info(f"Trajectory {r.__class__.__name__} start")
            r.on_trajectory_frame(self)
            self.info(f"Trajectory {r.__class__.__name__} finished")
        self.trajectory_frame += 1

    def get_openmm_topology(self) -> mmapp.Topology:
        """Helper that generates an OpenMM Topology object required for creating an OpenMMApp Simulation object.

        Chains will be set to initial molecules, residues will be set to initial residues.
        The periodic box, bonds, atom names and charges reflect the current state of the system.
        It's recommended to use other methods to query those though, as the purpose of this function is to
        facilitate using Martini Daemon as a .top file parser and then continue simulating in vanilla OpenMM.
        """
        top = mmapp.Topology()

        if self.__context is not None:
            _, box = self.__context.get_positions()

            top.setPeriodicBoxVectors([box.a, box.b, box.c])

        init_molecules = [
            (name, count, len(self.system.molecule_types[name].atoms))
            for name, count in self.system.initial_molecules
        ]

        # iteration over initial molecules
        # chains correspond to initial molecules
        # residues correspond to residues
        c_mol_type = 0
        c_mol_idx = 0
        c_index_in_mol = 0
        c_chain = top.addChain()
        last_residue = -1
        c_res = None
        must_be_new_residue = True
        atoms = []
        for i in range(self.system.num_atoms()):
            resid = self.system.get_res_id(i)
            if last_residue < resid:
                last_residue = resid
                res_name = self.system.get_res_name(i)
                c_res = top.addResidue(res_name, c_chain)
                must_be_new_residue = False
            else:
                assert not must_be_new_residue
            atoms.append(
                top.addAtom(
                    self.system.get_name(i),
                    None,
                    c_res,
                    formalCharge=self.system.get_charge(i),
                )
            )
            c_index_in_mol += 1
            if c_index_in_mol >= init_molecules[c_mol_type][2]:
                c_mol_idx += 1
                c_index_in_mol = 0
                c_chain = top.addChain()
                must_be_new_residue = True
            if c_mol_idx >= init_molecules[c_mol_type][1]:
                c_mol_idx = 0
                c_mol_type += 1

        for i, j in self.system.collect_bonds(
            ["bond", "constraint", "vsite"]
        ).to_list():
            top.addBond(atoms[i], atoms[j])

        return top

    def step(self, n_steps: int, traj=False, dm=False, silent=False):
        """Does the following:
        - steps n_steps
        - D/M algorithm if dm is true
        - trigger a trajectory frame on reporters if traj is True
        - display info to logs and screen, % info given by self.current_step and self.md_steps
        - update self.current_step

        Note: reinitializing is now handled automatically by context.

        Note: if calling step() manually, must manually call finish() after to properly flush output files!
        """
        start_time = time()
        if traj:
            self.__do_traj_frame()

        self.info(f"doing md steps to go from {self.current_step} to")
        self.current_step += n_steps
        percent = (
            self.current_step / self.total_steps * 100.0
            if self.total_steps > 0
            else 100.0
        )
        self.info(f"step {self.current_step}")
        if n_steps > 0:
            self.info(f"md_steps {n_steps}")
            self.info("Reinitialize start")
            self.context.do_steps(0)
            self.info("Reinitialize finished")
            self.info("MD start")
            self.context.do_steps(n_steps)
            self.info("MD finished")
        if self.total_steps > 0 and n_steps > 0 and not silent:
            self.time_ps += self.dt_ps * n_steps
            time_left = self.__format_time(
                self.__last_step_time * (self.total_steps - self.current_step)
            )
            reporter_data = " ".join(
                filter(None, [r.interactive_line(self) for r in self.__reporters])
            )
            sys.stdout.write(
                f"\033[2K\rstep {self.current_step}"
                f"({self.__format_sim_time(self.time_ps)}, "
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
                for r in self.__reporters:
                    r.pre_modification(self)
                # we now need copies of fragments, since they possibly got consumed in the reaction
                self.info("Modification start")
                reactions: list[tuple[str, list[Fragment]]] = self.top.modification(
                    reactions
                )
                self.info(f"After modification there were {len(reactions)} reactions")
                self.info("Modification finished")
                if len(reactions) > 0:
                    for r in self.__reporters:
                        self.info(f"OnReaction {r.__class__.__name__} start")
                        # local minimization happens here for example
                        r.on_reaction(self, reactions)
                        self.info(f"OnReaction {r.__class__.__name__} finished")
                    for r in self.__reporters:
                        self.info(f"PostReaction {r.__class__.__name__} start")
                        r.post_reaction(self)
                        self.info(f"PostReaction {r.__class__.__name__} finished")

        end_time = time()
        if n_steps > 0:
            step_time = (end_time - start_time) / n_steps
            if self.__last_step_time > 0.0:
                self.__last_step_time = step_time * 0.01 + self.__last_step_time * 0.99
            elif not traj or not dm:
                # step 0 tends to have both xtc and dm as True, and is
                # usually unrepresentatively slow
                scale = self.dm_frequency / self.traj_frequency
                if scale > 1.0:
                    scale = 1.0 / scale
                if self.__first_step_time == 0.0:
                    # continuations might not start with an expensive step
                    scale = 0.0
                self.__last_step_time = (
                    step_time * (1 - scale) + scale * self.__first_step_time
                )
            else:
                # step 0 probably
                self.__first_step_time = step_time
