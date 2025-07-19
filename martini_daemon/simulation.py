"""Daemon Simulation object (prototype)
the D/M algorithm + wrappers
"""

from .forces.nonbonded import NonBonded
from .top_parser import DaemonTopFile
from .meta import alias
from .gro_file import read_gro
from .utils import backup_try
from .reporters.checkpoint_reporter import load_checkpoint
import sys
import openmm as mm
import logging
from random import random
import time
from typing import Any
import math
import os


class Simulation():

    logger_id = 0

    @alias({
        "xtc_frequency": "traj_frequency",
        "gro_path": "geom_path"
    })
    def __init__(
        # main input files
        self, top_path: str = "", geom_path: str = "",
        # simulation length / frequencies
        md_steps: int = 0, dm_frequency: int = 0,
        traj_frequency: int = 5000,
        # other input / output files
        sim_name: str = "out",
        chk_path: str = "",
        restraint_coord_path=None,
        # integrator, context params, coupling
        integrator: mm.Integrator = None,
        coupling: list[mm.Force] = None,
        remove_com_motion: bool = True,
        platform: str | None | mm.Platform = None,
        context_parameters: None | dict[str, str] = None,
        nonbonded_force: NonBonded = None,
        # parsing things
        include_dir: str | None = None,
        defines: dict[str, str] = {},
        reporters: list[Any] = [],
        # T* additional params
        neighbor_cutoff: float = 1.1,
        force_reinitialize: bool = False,
        max_absolute_rate: float | None = None,
        rate_highest_probability: float = 1.0,
    ):
        """
        Input topology and geometry can be either:
            top_path -> .top file
            geom_path -> .gro file
            gro_path -> alias for geom_path
        or
            chk_path -> .chk file (simulation checkpoint)
        Both cannot be specified! Note: when loading from checkpoints
        include_dir and #define don't matter, they only matter for
        .top files.

        md_steps -> number of total steps to do
        dm_frequency -> how often to do a detection/modification
        traj_frequency -> how often to write trajectory
            (and associated reporters)
        xtc_frequency -> alias for traj_frequency
        sim_name -> prefix for simulation output files
        integrator -> integrator to use, by default LangevinMiddleIntegrator,
            0.02 ps timestep, 1 ps^-1 friction, 300K temp
        coupling -> additional openmm forces to add to the system.
        Useful for T/p coupling. By default contains a monte carlo barostat
        set to 300K temp and 1 bar.
        remove_com_motion -> should we remove center of mass motion?
        platform -> openmm platform, can be a string or obj
        minimize_energy -> should we minimize energy?
        generate_velocities -> should we generate velocities?
        include_dir -> search for .itp files here too
        defines -> #defines for .itp
        reporters -> list of Reporters
        neighbor_cutoff -> neighbor cutoff for the detection algorithm
        force_reinitialize -> only used to benchmark the impact of
            reinitialize, keep it False
        max_absolute_rate -> highest rate of reaction possible for rate limited
            reactions, for relative rate = 1, concentration = 1 per simulation box,
            in units of once per d/m frequency (units up to change)

            use case: slowing down reactions only at times where all reactions
            in the system are fast (relative to their concentration),
            without slowing them down in other cases
        highest_probability -> 0 to 1., for every rate controlled reaction,
            the probability of being accepted can be scaled by a value
            affects the general speed of reactions in most cases
        """

        # Self initialization
        use_checkpoint = False
        if geom_path == "" and top_path == "" and chk_path != "":
            use_checkpoint = True
        elif chk_path == "" and geom_path != "" and top_path != "":
            pass
        else:
            raise ValueError(
                "Please specify: geom_path and top_path only OR chk_path only."
            )
        if restraint_coord_path is None:
            restraint_coord_path = geom_path
        self.i: int = 0
        self.md_steps: int = md_steps
        self.traj_path: str = sim_name + ".xtc"
        self.out_path: str = sim_name + ".gro"
        self.log_path: str = sim_name + ".log"
        self.reactions: int = 0
        self.last_step_time: float = 0.
        self.first_step_time: float = 0.
        self.xtc_freq: int = traj_frequency
        self.dm_freq: int = dm_frequency
        self.neighbor_cutoff = neighbor_cutoff
        if type(platform) is str:
            platform = mm.Platform.getPlatformByName(platform)
        include_dir = include_dir or (
            "GMXDATA" in os.environ and
            os.path.join(os.environ["GMXDATA"], "top")
        ) or (
            "GMXBIN" in os.environ and
            os.path.join(os.environ["GMXBIN"], "..", "share", "gromacs", "top")
        ) or "/usr/local/gromacs/share/gromacs/top"
        self.force_reinitialize = force_reinitialize
        if integrator is None:
            integrator = mm.LangevinMiddleIntegrator(
                300 * mm.unit.kelvin,
                1.0 / mm.unit.picosecond,
                0.02 * mm.unit.picosecond
            )
        if coupling is None:
            coupling = [
                mm.MonteCarloBarostat(
                    1.0 * mm.unit.bar,
                    300 * mm.unit.kelvin
                )
            ]
        if nonbonded_force is None:
            nonbonded_force = NonBonded(epsilon_r=15.0, cutoff_nm=1.1)
        self.dt_ns: float = (
            integrator.getStepSize().value_in_unit(mm.unit.nanosecond)
        )

        # Logging setup
        backup_try(self.log_path)
        self.logger = logging.getLogger(f"logger_{self.logger_id}_{random()}")
        self.logger_id += 1
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s %(message)s")
        fh = logging.FileHandler(self.log_path)
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)
        stream_formatter = logging.Formatter("[%(levelname)s] %(message)s")
        sh = logging.StreamHandler(sys.stderr)
        sh.setFormatter(stream_formatter)
        sh.setLevel(logging.WARN)
        self.logger.addHandler(sh)
        self.logger.info("DaemonSimulation __init__ called")
        self.logger.info(f"Parameters: {top_path} {geom_path} "
                         f"dt (ns): {self.dt_ns} "
                         f"steps: {md_steps} dm_freq: {dm_frequency} "
                         f"traj_freq: {traj_frequency} "
                         f"sim_name: {sim_name} platform: {platform}")

        self.logger.info("Parsing start")
        # Parsing - No checkpoint
        if not use_checkpoint:
            _, respos, _ = read_gro(restraint_coord_path)
            self.system, self.top = DaemonTopFile(
                top_path, nonbonded_force,
                include_dir=include_dir, defines=defines,
                nlist_cutoff=neighbor_cutoff,
                max_absolute_rate=max_absolute_rate,
                rate_highest_probability=rate_highest_probability,
                logger=self.logger,
                respos=respos
            )
            # benefits of function based scope
            box, pos, vel = read_gro(geom_path)
        # Parsing - Checkpoint
        else:
            if (
                restraint_coord_path is not None and \
                len(restraint_coord_path) > 0
            ):
                raise ValueError(
                    "Currently position restraints don't work from checkpoints"
                )
            self.i, self.system, self.top = load_checkpoint(
                chk_path, self.logger, nonbonded_force,
                neighbor_cutoff,
                max_absolute_rate, rate_highest_probability
            )
        self.logger.info("Parsing finished")

        # Coupling and integrators
        self.logger.info("Setup integrator and coupling start")
        for f in coupling:
            self.system.add_force(f)
        if remove_com_motion:
            self.system.remove_com_motion()

        # Context build and reporter initialization
        self.logger.info("Context build start")
        if use_checkpoint:
            self.system.load_finish(
                integrator, platform, context_parameters
            )
        else:
            self.system.build_context(
                integrator, box, platform, context_parameters
            )
            self.system.set_positions(pos)
            if vel is not None:
                self.system.set_velocities(vel)

        for rep in reporters:
            rep._sysstar = self.system
            rep._topstar = self.top
            rep._simulation = self
            self.system.add_reporter(rep)
            self.top.add_reporter(rep)
        self.reporters = reporters

        self.system.set_xtc_path(self.traj_path)
        self.top.init_dm(sim_name)
        if not load_checkpoint:
            self.system.write_xtc_frame(0, 0)
            self.logger.info("Initial geometry XTC frame written")
        self.logger.info("Context build finished")
        self.logger.info("__init__ end")

    def generate_velocities(self, T=300):
        """Generate velocities at temp T (in kelvin)."""
        self.logger.info("genvel start")
        self.system.generate_velocities(T)
        self.logger.info("genvel finished")

    def minimize_energy(
        self, tolerance=10, max_steps=0, write_xtc=False, out=None
    ):
        """
        Minimizes energy of the system. If write_xtc is True, an xtc frame
        will be written. If out is set to a string value, a .gro file of
        the minimized coordinates will be written there.
        """
        self.logger.info("Energy min start")
        self.system.minimize_energy(tolerance, max_steps)
        if write_xtc:
            self.system.write_xtc_frame(0, 0)
            self.logger.info("Energy minimized XTC frame written")
        if out is not None:
            self.system.write_gro(out)
        self.logger.info("Energy min finished")

    def simulate(self):
        remaining = self.md_steps - self.i
        sim_ns = self.md_steps * self.dt_ns
        print(f"Simulation of {remaining} steps ({sim_ns} ns)")
        gcd = math.gcd(self.xtc_freq, self.dm_freq, remaining)
        print(f"D/M freq {self.dm_freq} XTC freq {self.xtc_freq} gcd {gcd}")
        while self.i < self.md_steps:
            self.step(
                gcd, xtc=self.i % self.xtc_freq == 0 if self.xtc_freq > 0 else False,
                dm=self.i % self.dm_freq == 0 if self.dm_freq > 0 else False
            )
        print()
        self.system.write_gro(self.out_path)
        for reporter in self.reporters:
            reporter.finish()

    def step(self, steps=1, xtc=True, dm=True):
        """
            Do a step of the following:
            - steps MD steps
            - D/M algorithm if dm is true
            - reinitialize system if reactions happened
            - write an XTC frame if xtc is True
            - display info to logs and screen, % info given by self.i and self.md_steps
            - update self.i
        """
        start_time = time.time()
        self.i += steps
        percent = self.i/self.md_steps*100 if self.md_steps > 0 else 100
        self.logger.info(f"step {self.i}")
        if steps > 0:
            self.logger.info(f"md_steps {steps}")
            self.logger.info("MD start")
            self.system.do_steps(steps)
            self.logger.info("MD finished")
        if self.md_steps > 0 and steps > 0:
            ns_so_far = self.dt_ns * self.i
            time_left = self.last_step_time * (self.md_steps - self.i)
            time_fmt: str
            if time_left < 3600:
                time_fmt = time.strftime("%M:%S", time.gmtime(time_left))
            elif time_left < 3600 * 24:
                time_fmt = time.strftime("%H:%M:%S", time.gmtime(time_left))
            elif time_left < 3600 * 24 * 30:
                time_fmt = time.strftime("%dd %H:%M:%S", time.gmtime(time_left))
            else:
                time_fmt = f"Longer than a month ({time_left} seconds)"
            reporter_data = " ".join(filter(None, [r.interactive_line() for r in self.reporters]))
            sys.stdout.write(f"\033[2K\rstep {self.i}"
                             f"({ns_so_far:.2f} ns, "
                             f"{percent:.1f}%) "
                             f"{time_fmt} {reporter_data}")
        if dm:
            self.logger.info("Detection start")
            pos, box = self.system.get_positions()
            reactions = self.top.detection(self.i, box, pos)
            self.logger.info("Detection finished")
            if len(reactions) > 0:
                self.logger.info("Modification start")
                self.reactions += len(reactions)
                self.top.modification(self.i, reactions)
                self.logger.info("Modification finished")
            if len(reactions) > 0 or self.force_reinitialize:
                self.logger.info("Reinitialize start")
                self.system.reinitialize(self.force_reinitialize)
                self.logger.info("Reinitialize finished")
            self.logger.info(f"reactions {len(reactions)}")

        if xtc:
            self.logger.info("XTC write start")
            self.system.write_xtc_frame(self.i, self.xtc_freq)
            self.logger.info("XTC write finished")
        end_time = time.time()
        # timing info update
        if steps > 0:
            step_time = (end_time - start_time) / steps
            if self.last_step_time > 0.:
                self.last_step_time = step_time * 0.01 + self.last_step_time * 0.99
            elif not xtc or not dm:
                # step 0 tends to have both xtc and dm as True, and is
                # usually unrepresentatively slow
                scale = self.dm_freq / self.xtc_freq
                if scale > 1.:
                    scale = 1. / scale
                if self.first_step_time == 0.:
                    # continuations might not start with an expensive step
                    scale = 0.
                self.last_step_time = step_time * (1 - scale) + scale * self.first_step_time
            else:
                # step 0 probably
                self.first_step_time = step_time
