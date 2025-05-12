"""Daemon Simulation object (prototype)
the D/M algorithm + wrappers
"""

from .top_parser import DaemonTopFile
from .sysstar import SysStar
from .topstar import TopStar
from .meta import alias
from .gro_file import read_gro
from .utils import backup_try, smooth
from .reporters.checkpoint_reporter import load_checkpoint
from .reporters.reporter import Reporter
import sys
import openmm as mm
import logging
from random import random
import time
from typing import Any
import math
import os


class DaemonSimulation():

    system: SysStar
    top: TopStar

    i: int
    logger: logging.Logger
    logger_id = 0
    reactions: int
    last_step_time: tuple[float, float]
    md_steps: int
    dm_freq: int
    xtc_freq: int
    traj_path: str
    out_path: str
    log_path: str
    neighbor_cutoff: float
    dt_ns: float
    force_reinitialize: bool
    reporters: list[Reporter]
    smooth_params = (0.002, 0.01)

    @alias({
        "xtc_frequency": "traj_frequency",
        "gro_path": "geom_path",
        "T_kelvin": "t_kelvin",
        "T_type": "t_type"
    })
    def __init__(self, top_path: str = "", geom_path: str = "",
                 md_steps: int = 0, dm_frequency: int = 0,
                 chk_path: str = "",
                 traj_frequency: int = 5000,
                 sim_name: str = "out",
                 t_kelvin: float | None = 300., t_type: str = "langevin",
                 p_bar: float = 1.,
                 dt_ps: float = 0.02,
                 platform: str | None | mm.Platform = None,
                 minimize_energy: bool = True,
                 generate_velocities: bool = True,
                 remove_com_motion: bool = True,
                 epsilon_r: float = 15.0,
                 nonbonded_cutoff_nm: float = 1.1,
                 include_dir: str | None = None,
                 defines: dict[str, str] = {},
                 reporters: list[Any] = [],
                 friction_ps_1: float = 2.0,
                 neighbor_cutoff: float = 1.1,
                 force_reinitialize: bool = False,
                 max_absolute_rate: float | None = None,
                 rate_smoothing: tuple[float, float] = (0.01, 0.02),
                 rate_highest_probability: float = 1.0
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
        t_kelvin -> reference temperature
        T_kelvin -> alias for t_kelvin
        t_type -> "langevin" or "andersen"
        T_type -> alias for t_type
        p_bar -> reference pressure in bars
        dt_ps -> single MD timestep in ps
        platform -> openmm platform
        minimize_energy -> should we minimize energy?
        generate_velocities -> should we generate velocities?
        remove_com_motion -> should we remove center of mass motion?
        epsilon_r
        nonbonded_cutoff_nm -> LJ / ES neighbor list cutoff
        include_dir -> search for .itp files here too
        defines -> #defines for .itp
        reporters -> list of Reporters
        friction_ps_1 -> friction value for temperature coupling
        neighbor_cutoff -> neighbor cutoff for the detection algorithm
        force_reinitialize -> only used to benchmark the impact of
            reinitialize, keep it False
        max_absolute_rate -> highest rate of reaction possible for rate limited
            reactions, for relative rate = 1, concentration = 1 per simulation box,
            in units of once per d/m frequency (units up to change)

            use case: slowing down reactions only at times where all reactions
            in the system are fast (relative to their concentration),
            without slowing them down in other cases
        rate_smoothing_constant -> parameters for the smoothing algorithm,
            the smaller the slower it is to react to rate changes. The first
            argument is for value smoothing, the second is for slope smoothing
            (it uses double exponential smoothing)
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
        self.i: int = 0
        self.md_steps: int = md_steps
        self.traj_path: str = sim_name + ".xtc"
        self.out_path: str = sim_name + ".gro"
        self.log_path: str = sim_name + ".log"
        self.reactions: int = 0
        self.last_step_time: tuple[float, float] = (0., 0.)
        self.first_step_time: float = 0.
        self.xtc_freq: int = traj_frequency
        self.dm_freq: int = dm_frequency
        self.neighbor_cutoff = neighbor_cutoff
        if t_kelvin is None and t_type != "none":
            raise ValueError("No valid T temperature given")
        if t_type not in {"andersen", "langevin", "none"}:
            raise ValueError("Unknown t_type")
        T = t_kelvin * mm.unit.kelvin if t_type != "none" else None
        p = p_bar * mm.unit.bar if p_bar is not None else None
        if p is None and t_type == "none":
            raise ValueError("Must couple T for p coupling")
        dt = dt_ps * mm.unit.picosecond
        self.dt_ns: float = dt.value_in_unit(mm.unit.nanosecond)
        nonbonded_cutoff = nonbonded_cutoff_nm * mm.unit.nanometer
        friction = friction_ps_1 / mm.unit.picosecond
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
                         f"T: {T} p: {p} dt: {dt} "
                         f"steps: {md_steps} dm_freq: {dm_frequency} "
                         f"traj_freq: {traj_frequency} "
                         f"sim_name: {sim_name} platform: {platform}")

        self.logger.info("Parsing start")
        # Parsing - No checkpoint
        if not use_checkpoint:
            self.system, self.top = DaemonTopFile(
                top_path,
                include_dir=include_dir, defines=defines,
                epsilon_r=epsilon_r, nonbonded_cutoff=nonbonded_cutoff,
                nlist_cutoff=neighbor_cutoff,
                max_absolute_rate=max_absolute_rate,
                rate_smoothing=rate_smoothing,
                rate_highest_probability=rate_highest_probability,
                logger=self.logger
            )
            # benefits of function based scope
            box, pos, vel = read_gro(geom_path)
        # Parsing - Checkpoint
        else:
            self.i, self.system, self.top = load_checkpoint(
                chk_path, self.logger, epsilon_r, nonbonded_cutoff,
                neighbor_cutoff
            )
        self.logger.info("Parsing finished")

        # Coupling and integrators
        self.logger.info("Setup integrator and coupling start")
        if t_type == "andersen":
            self.system.add_force(mm.AndersenThermostat(T, friction))

        if p is not None:
            self.system.add_force(mm.MonteCarloBarostat(p, T))
        if remove_com_motion:
            self.system.add_force(mm.CMMotionRemover())

        integrator: mm.Integrator
        if t_type == "langevin":
            integrator = mm.LangevinIntegrator(T, friction, dt)
        else:
            integrator = mm.VerletIntegrator(dt)
        self.logger.info("Setup integrator and coupling finished")

        # Context build and reporter initialization
        self.logger.info("Context build start")
        if use_checkpoint:
            self.system.load_finish(integrator, platform)
        else:
            self.system.build_context(integrator, box, platform)
            self.system.set_positions(pos)

            if not generate_velocities:
                if vel is None:
                    raise ValueError("Generate velocities is false, but there are no velocities in the .gro file.")
                if minimize_energy:
                    raise ValueError("Attempt to minimize energy and then read velocities from a file, which is likely wrong.")
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

        # Genvel, energy min
        if use_checkpoint and (generate_velocities or minimize_energy):
            raise ValueError("Cannot load .chk and generate velocities or minimize energies.")
        if generate_velocities:
            self.logger.info("genvel start")
            self.system.generate_velocities(T)
            self.logger.info("genvel finished")
        if minimize_energy:
            self.logger.info("Energy min start")
            self.system.minimize_energy()
            self.system.write_xtc_frame(0, 0)
            self.logger.info("Energy minimized XTC frame written")
            self.logger.info("Energy min finished")

        self.logger.info("__init__ end")

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

    def step(self, steps=1, xtc=True, dm=True, neighbor=True):
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
        self.logger.info(f"md_steps {steps}")
        self.logger.info("MD start")
        self.system.do_steps(steps)
        self.logger.info("MD finished")
        if self.md_steps > 0:
            ns_so_far = self.dt_ns * self.i
            time_left = self.last_step_time[0] * (self.md_steps - self.i)
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
        step_time = (end_time - start_time) / steps
        if self.last_step_time[0] > 0.:
            self.last_step_time = smooth(step_time, self.last_step_time, self.smooth_params)
        elif not xtc or not dm:
            # step 0 tends to have both xtc and dm as True, and is
            # usually unrepresentatively slow
            scale = self.dm_freq / self.xtc_freq
            if scale > 1.:
                scale = 1. / scale
            if self.first_step_time == 0.:
                # continuations might not start with an expensive step
                scale = 0.
            self.last_step_time = (step_time * (1 - scale) + scale * self.first_step_time, 0.)
        else:
            # step 0 probably
            self.first_step_time = step_time
