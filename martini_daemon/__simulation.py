"""Friendly Simulation API."""

import math
import os
import shutil
import sys
import traceback
from collections.abc import Callable
from datetime import datetime
from importlib.metadata import version
from time import time
from typing import Any

import numpy as np
import numpy.typing as npt
import openmm as mm
import openmm.app as mmapp
from openmm.unit import md_unit_system

from .__core import Context, System, wrap_coupling
from .__forces import NonBonded
from .__formats import Checkpoint, read_geometry, write_geometry
from .__parser import GromacsTopFile, InvalidTopologyError
from .__reporter import Reporter
from .__rust import Fragment, PeriodicBox, build_version
from .__topstar import TopStar


def _format_time(total_seconds: float | int) -> str:
    c = int(total_seconds)
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


def _format_sim_time(ps: float) -> str:
    if ps < 1000.0:
        return f"{ps:.2f} ps"
    if ps < 1000000.0:
        return f"{ps / 1000.0:.2f} ns"
    return f"{ps / 1000000:.2f} μs"


class Simulation:
    """Simulation class."""

    def __init__(
        self,
        topology: str,
        geometry: str
        | None
        | tuple[PeriodicBox, npt.NDArray[np.float64], npt.NDArray[np.float64] | None]
        | Checkpoint,
        md_steps: int,
        reporters: list[Reporter] | None = None,
        dm_frequency: int = 0,
        traj_frequency: int = 0,
        sim_name: str = "out",
        coupling: list[mm.Force] | None = None,
        integrator: mm.Integrator | None = None,
        options: dict[str, Any] | None = None,
        include_dirs: list[str] | None = None,
        defines: dict[str, str] | None = None,
        platform: str | None | mm.Platform = None,
        context_parameters: None | dict[str, str] = None,
        nonbonded: Callable[[System], NonBonded] | type[NonBonded] | None = None,
        continue_sim: bool = False,
    ) -> None:
        """Create a simulation.

        * Provides a friendly interface for reporters requesting output files. Contains default values for Martini simulations.
        * Holds simulation metadata, such as current step, simulation name.
        * Owns all reporters and the log file handle.
        * Is passed around to all reporters to provide the required metadata for reporting.
        * Is the required glue between all components, also only uses the public interface of different components.
        * Provides access to the system, context, topstar instance to all reporters and the user.

        :param topology: Path to the Martini .top file.
        :param geometry: Path to the geometry file (.gro, .xyz). Alternatively, a tuple of PeriodicBox, positions
            (numpy array) and velocities (numpy array, optional, can be None) is also accepted.
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
        self.total_steps: int = md_steps
        self.__sim_name = sim_name
        # metadata
        if isinstance(geometry, Checkpoint):
            self.trajectory_frame = geometry.trajectory_frame
            self.reactions_so_far = geometry.reactions_so_far
            self.current_step = geometry.current_step
            self.time_ps = geometry.time_ps
        else:
            self.trajectory_frame: int = 0
            self.reactions_so_far: int = 0
            self.current_step: int = 0
            self.time_ps: float = 0.0

        # to avoid double finish
        self.__finished = False

        # for the estimated time left display
        self.__last_step_time = 0.0
        self.__first_step_time = 0.0

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

        self.log = open(self.request_path(".log", copy=continue_sim), "a")
        self.info(f"Martini Daemon {version('martini_daemon')} log file")
        self.info("Build version:", build_version())
        self.info(
            "Parameters:",
            topology,
            geometry,
            f"steps: {self.total_steps}",
            f"dm_freq: {self.dm_frequency}",
            f"traj_freq: {self.traj_frequency}",
            f"sim_name: {self.__sim_name}",
            f"platform: {platform}",
            f"context_parameters: {context_parameters}",
            f"defines: {defines}",
            f"include_dirs: {include_dirs}",
        )
        if options is None:
            options = {}
        if options.get("epsilon_r") is None:
            options["epsilon_r"] = 15.0
        if options.get("cutoff") is None:
            options["cutoff"] = 1.1
        if options.get("respos") is None and type(geometry) is str:
            self.info(f"Loading restraint reference positions from {geometry}.")
            options["respos"] = read_geometry(geometry)[1]

        # Parsing
        self.info("Parsing start")
        self.system: System = System(options=options)
        try:
            GromacsTopFile(
                self.system, topology, include_dirs=include_dirs, defines=defines
            )
        except InvalidTopologyError:
            self.error("Fatal error during .top parsing.")
            # don't print full traceback
            raise InvalidTopologyError from None
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
            f"dt (ps): {self.dt_ps}",
            f"type: {type(md_integrator).__name__}",
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
        if geometry is not None:
            if type(geometry) is str:
                box, start_pos, start_vel = read_geometry(geometry)
                self.info(
                    f"Read {len(start_pos)} atoms from {geometry}. Box: {box.to_lattice()}. Velocities read? {start_vel is not None}."
                )
            elif isinstance(geometry, Checkpoint):
                box = geometry.box
                start_pos = geometry.pos
                start_vel = geometry.vel
            else:
                box, start_pos, start_vel = geometry
                self.info("Initial box, pos, vel provided as a tuple.")

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
            Generally, only pass True, if you intend to append to the contents.
        """
        path = self.__sim_name + suffix
        self.__backup_try(path, copy=copy)
        return path

    def finish(self) -> None:
        """End the simulation.

        Runs finish on all reporters. Reporters should close their own file handles.

        Called by simulate(), or should be called by the user manually otherwise.
        """
        if self.__finished:
            return
        self.__finished = True
        for r in self.__reporters:
            r.on_simulation_finish(self)
        self.log.close()

    def info(self, *args: str) -> None:
        """Write a message to log."""
        self.log.write(
            " ".join(
                [datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S,%f")[:-3], *args]
            )
        )

    def warn(self, message: str) -> None:
        """Write a warning to the log file and to stderr."""
        self.info("[WARNING] " + message)
        print("[WARNING]", message, file=sys.stderr)

    def error(self, message: str) -> None:
        """Write an error to the log file and stderr."""
        self.info("[ERROR] " + message)
        print("[ERROR]", message, file=sys.stderr)

    # Friendly interface for setting up and running simulations
    @classmethod
    def set_process_title(cls, newname: bytes = b"daemon") -> None:
        """Set process title to something else than "python".

        Nothing critical, purely aesthetic.
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

    def save_geometry(self, path: str) -> None:
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

    def simulate(self, finish: bool = True) -> None:
        """Run the simulation.

        Performs the remaining steps (self.total_steps - self.current_step),
        and then calls self.finish() if finish is True.

        Will perform Detection/Modification and Trajectory writing according to their frequencies specified in
        the constructor for Simulation. Catches uncaught exceptions and logs them,
        and finishes the simulation prematurely if one occurs.

        :param finish: whether to close all output files. If False, must call finish() manually.
        """
        remaining = self.total_steps - self.current_step
        sim_ps = self.total_steps * self.dt_ps
        print(f"Simulation of {remaining} steps ({_format_sim_time(sim_ps)})")
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

    def __do_traj_frame(self) -> None:
        """Write a frame to all trajectory files."""
        for r in self.__reporters:
            self.info(f"Trajectory {r.__class__.__name__} start")
            r.on_trajectory_frame(self)
            self.info(f"Trajectory {r.__class__.__name__} finished")
        self.trajectory_frame += 1

    def get_openmm_topology(self) -> mmapp.Topology:
        """Get an OpenMM Topology object.

        Helper that generates an OpenMM Topology object required for creating an OpenMMApp Simulation object.

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

    def step(
        self, n_steps: int, traj: bool = False, dm: bool = False, silent: bool = False
    ) -> None:
        """Perform simulation steps.

        Does the following:
        - Write a trajectory frame
        - The integrator does n_steps MD steps
        - Update self.current_step and self.time_ps
        - Process reactions based on the Detection/Modification (D/M) algorithm
        - Display info to logs and screen, % info given by self.current_step and self.md_steps

        :param n_steps: number of MD steps
        :param traj: whether to write a trajectory frame
        :param dm: whether to run the D/M algorithm
        :param silent: whether to suppress printing a status line to stdout

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
            time_left = _format_time(
                self.__last_step_time * (self.total_steps - self.current_step)
            )
            reporter_data = " ".join(
                filter(None, [r.interactive_line(self) for r in self.__reporters])
            )
            sys.stdout.write(
                f"\033[2K\rstep {self.current_step}"
                f"({_format_sim_time(self.time_ps)}, "
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
                    self.reactions_so_far += len(reactions)

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
