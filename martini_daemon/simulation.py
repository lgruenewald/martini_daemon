#!/usr/bin/env python3
"""Daemon Simulation object (prototype)
the D/M algorithm + wrappers
"""

from .top_parser import DaemonTopFile
from .sysstar import SysStar
from .topstar import TopStar, ReactionTemplate, Fragment
from .reporters.reporter import Reporter
import sys
import openmm as mm  # type: ignore[import-untyped]
from openmm.app import GromacsGroFile  # type: ignore[import-untyped]
from .utils import backup_try
import logging
from random import random
import time
from typing import Any
import math
import os


class DaemonSimulation():

    system: SysStar
    top: TopStar
    gro: GromacsGroFile

    reaction_list: list[ReactionTemplate]
    init_list: list[Fragment]
    init_map: dict[str, list[int]]

    logger_id = 0
    reactions: int
    md_steps: int
    dm_freq: int
    xtc_freq: int
    traj_path: str  # trajectory to write
    out_path: str  # final geometry to write

    def __init__(self, top_path: str, gro_path: str,
                 md_steps: int, dm_frequency: int,
                 xtc_frequency: int = 5000,
                 sim_name: str = "out",
                 T_kelvin: float | None = 300., T_type: str = "langevin",
                 p_bar: float = 1., dt_ps: float = 0.02,
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
                 neighbor_cutoff: float = 1.5,
                 ):

        # Self initialization
        self.md_steps: int = md_steps
        self.traj_path: str = sim_name + ".xtc"
        self.out_path: str = sim_name + ".gro"
        self.log_path: str = sim_name + ".log"
        self.reactions: int = 0
        self.last_step_time: int = 0.
        self.xtc_freq: int = xtc_frequency
        self.dm_freq: int = dm_frequency
        self.neighbor_cutoff = neighbor_cutoff
        if T_kelvin is None and T_type != "none":
            raise ValueError("No valid T temperature given")
        if T_type not in {"andersen", "langevin", "none"}:
            raise ValueError("Unknown T_type")
        T = T_kelvin * mm.unit.kelvin if T_type != "none" else None
        p = p_bar * mm.unit.bar if p_bar is not None else None
        if p is None and T_type == "none":
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
        self.logger.info(f"Parameters: {top_path} {gro_path} "
                         f"T: {T} p: {p} dt: {dt} "
                         f"steps: {md_steps} dm_freq: {dm_frequency} "
                         f"xtc_freq: {xtc_frequency} "
                         f"sim_name: {sim_name} platform: {platform}")

        # Parsing
        self.logger.info("Parsing start")
        self.system, self.top = DaemonTopFile(
            top_path,
            include_dir=include_dir, defines=defines,
            epsilon_r=epsilon_r, nonbonded_cutoff=nonbonded_cutoff,
            nlist_cutoff=neighbor_cutoff,
            logger=self.logger
        )
        self.gro = GromacsGroFile(gro_path)
        box = self.gro.getPeriodicBoxVectors()
        pos = self.gro.getPositions(asNumpy=True)
        self.logger.info("Parsing finished")

        # Coupling and integrators
        self.logger.info("Setup integrator and coupling start")
        if T_type == "andersen":
            self.system.add_force(mm.AndersenThermostat(T, friction))

        if p is not None:
            self.system.add_force(mm.MonteCarloBarostat(p, T))
        if remove_com_motion:
            self.system.add_force(mm.CMMotionRemover())

        integrator: mm.Integrator
        if T_type == "langevin":
            integrator = mm.LangevinIntegrator(T, friction, dt)
        else:
            integrator = mm.VerletIntegrator(dt)
        self.logger.info("Setup integrator and coupling finished")

        # Context build and reporter initialization
        self.logger.info("Context build start")
        if platform is not None:
            self.system.build_context(integrator, box, platform)
        else:
            self.system.build_context(integrator, box)
        self.system.set_positions(pos)
#        self.system._context.setVelocities() TODO from the .gro file

        for rep in reporters:
            rep._sysstar = self.system
            rep._topstar = self.top
            self.system.add_reporter(rep)
            self.top.add_reporter(rep)

        self.system.set_xtc_path(self.traj_path)
        self.system.write_xtc_frame(0)
        self.logger.info("Initial geometry XTC frame written")
        self.logger.info("Context build finished")

        # Genvel, energy min
        if generate_velocities:
            self.logger.info("genvel start")
            self.system.generate_velocities(T)
            self.logger.info("genvel finished")
        if minimize_energy:
            self.logger.info("Energy min start")
            self.system.minimize_energy()
            self.system.write_xtc_frame(0)
            self.logger.info("Energy minimized XTC frame written")
            self.logger.info("Energy min finished")

        self.logger.info("__init__ end")

    def simulate(self):
        # Greatest common divisor
        sim_ns = self.md_steps * self.dt_ns
        print(f"Simulation of {self.md_steps} steps ({sim_ns} ns)")
        gcd = math.gcd(self.xtc_freq, self.dm_freq)
        print(f"D/M freq {self.dm_freq} XTC freq {self.xtc_freq} gcd {gcd}")
        for i in range(0, self.md_steps, gcd):
            self.step(
                gcd, xtc=i % self.xtc_freq == 0,
                dm=i % self.dm_freq == 0,
                i=i, max_steps=self.md_steps
            )
        print()
        self.system.write_gro(self.out_path)

    def step(self, steps=1, xtc=True, dm=True, neighbor=True, i=0, max_steps=0):
        """
            Do a step of the following:
            - steps MD steps
            - D/M algorithm if dm is true
            - reinitialize system if reactions happened
            - write an XTC frame if xtc is True
            - display info to logs and screen, % info given by i and max_steps
        """
        start_time = time.time()
        percent = i/max_steps*100 if max_steps > 0 else 100
        self.logger.info(f"step {i}")
        self.logger.info("MD start")
        self.system.do_steps(steps)
        self.logger.info("MD finished")
        if max_steps > 0:
            ns_so_far = self.dt_ns * i
            time_left = self.last_step_time * (max_steps - i)
            time_fmt: str
            if time_left < 3600:
                time_fmt = time.strftime("%M:%S", time.gmtime(time_left))
            elif time_left < 3600 * 24:
                time_fmt = time.strftime("%H:%M:%S", time.gmtime(time_left))
            elif time_left < 3600 * 24 * 30:
                time_fmt = time.strftime("%dd %H:%M:%S", time.gmtime(time_left))
            else:
                time_fmt = f"Longer than a month ({time_left} seconds)"
            sys.stdout.write(f"\033[2K\rstep {i}"
                             f"({ns_so_far:.2f} ns, "
                             f"{percent:.1f}%)\t"
                             f"{time_fmt}\t"
                             f"{self.reactions} reactions")
        if dm:
            self.logger.info("Detection start")
            pos, box = self.system.get_positions()
            reactions = self.top.dm_detection(
                i, box, pos
            )
            self.logger.info("Detection finished")
            if len(reactions) > 0:
                self.logger.info("Modification start")
                self.reactions += len(reactions)
                self.top.dm_modification(i, reactions)
                self.logger.info("Modification finished")
                self.logger.info("Reinitialize start")
                self.system.reinitialize()
                self.logger.info("Reinitialize finished")
        if xtc:
            self.logger.info("XTC write start")
            self.system.write_xtc_frame(self.xtc_freq)
            self.logger.info("XTC write finished")
        end_time = time.time()
        step_time = (end_time - start_time) / steps
        self.last_step_time = (
            0.99 * self.last_step_time + 0.01 * (step_time)
            if self.last_step_time > 0. else (step_time)
        )


if __name__ == "__main__":
    argv = sys.argv
    argc = len(sys.argv)

    if argc != 3:
        print("Usage: ./main.py <top file> <gro file>")
        quit(1)

    top_path = argv[1]
    gro_path = argv[2]

    sim = DaemonSimulation(top_path, gro_path)
    sim.simulate()
