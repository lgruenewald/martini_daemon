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
from openmm.unit import femtosecond, nanometer, nanosecond  # type: ignore[import-untyped]
from .utils import backup_try
import logging
import random
import time


class DaemonSimulation():

    system: SysStar
    top: TopStar
    gro: GromacsGroFile

    reaction_list: list[ReactionTemplate]
    init_list: list[Fragment]
    init_map: dict[str, list[int]]

    logger_id = 0
    reactions: int
    max_steps: int
    steps_per_step: int
    traj_path: str  # trajectory to write
    out_path: str  # final geometry to write

    def __init__(self, top_path, gro_path, T=300., p=1., dt=20*femtosecond,
                 max_steps=100, steps_per_step=5000,
                 sim_name="out", traj_path=None,
                 out_path=None, platform=None,
                 minimize_energy=True, generate_velocities=True,
                 remove_com_motion=True, epsilon_r=15.0,
                 nonbonded_cutoff=1.1*nanometer, include_dir=None,
                 defines={}, log_path=None, reporters=[],
                 friction=2.0, T_type="langevin", xtc_every=1):
        if traj_path is None:
            traj_path = sim_name + ".xtc"
        if out_path is None:
            out_path = sim_name + ".gro"
        if log_path is None:
            log_path = sim_name + ".log"
        backup_try(log_path)
        self.logger = logging.getLogger(f"logger_{self.logger_id}_{random.random()}")
        self.logger_id += 1
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s %(message)s")
        fh = logging.FileHandler(log_path)
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
                         f"max_steps: {max_steps} per_step {steps_per_step} "
                         f"sim_name: {sim_name} platform {platform}")
        self.logger.info("Parsing start")
        self.system, self.top = DaemonTopFile(
            top_path,
            include_dir=include_dir, defines=defines,
            epsilon_r=epsilon_r, nonbonded_cutoff=nonbonded_cutoff,
            logger=self.logger
        )
        self.logger.info("Parsing finished")
        self.logger.info("Coord read start")
        self.gro = GromacsGroFile(gro_path)
        self.logger.info("Coord read finished")
        if platform is not None:
            platform = mm.Platform.getPlatformByName(platform)

        if T is not None and T_type not in {"andersen", "langevin"}:
            raise ValueError("Unknown T_type")

        if T is not None and T_type == "andersen":
            self.system.add_force(mm.AndersenThermostat(T, friction))

        if p is not None:
            self.system.add_force(mm.MonteCarloBarostat(p, T))
        if remove_com_motion:
            self.system.add_force(mm.CMMotionRemover())

        integrator = None
        if T_type == "langevin":
            integrator = mm.LangevinIntegrator(T, friction, dt)
        else:
            integrator = mm.VerletIntegrator(dt)
        box = self.gro.getPeriodicBoxVectors()

        self.logger.info("Context build start")
        if platform is not None:
            self.system.build_context(integrator, box, platform)
        else:
            self.system.build_context(integrator, box)
        self.logger.info("Context build finished")

        self.logger.info("Setting positions")
        self.system.set_positions(self.gro.getPositions(True))
        for rep in reporters:
            # TODO only take instances
            if isinstance(rep, Reporter):
                rep._sysstar = self.system
                rep._topstar = self.top
            else:
                rep = rep(self.system, self.top)
            self.system.add_reporter(rep)
            self.top.add_reporter(rep)
        # must set xtc path after adding reporters currently
        self.logger.info("Writing initial positions to XTC")
        self.system.set_xtc_path(traj_path)
        self.system.write_xtc_frame(0)
        if generate_velocities:
            self.logger.info("genvel start")
            self.system.generate_velocities(T)
            self.logger.info("genvel finished")
        if minimize_energy:
            self.logger.info("Energy min start")
            self.system.minimize_energy()
            self.logger.info("Energy min finished")
            self.logger.info("Writing energy minimized positions to XTC")
            self.system.write_xtc_frame(0)
        self.max_steps = max_steps
        self.steps_per_step = steps_per_step
        self.traj_path = traj_path
        self.out_path = out_path
        self.reactions = 0
        self.last_step_time = 0.
        self.xtc_every = xtc_every
        self.step_ns = steps_per_step * dt.value_in_unit(nanosecond)
        self.logger.info("__init__ end")

    def simulate(self):
        for i in range(self.max_steps):
            self.step(i+1, self.max_steps)
        print()
        self.system.write_gro(self.out_path)

    def step(self, i=0, max_steps=0):
        """
            Do a step of the following:
            - self.steps_per_step MD steps
            - D/M algorithm
            - reinitialize system
            - write an XTC frame
        """
        start_time = time.time()
        percent = i/max_steps*100 if max_steps > 0 else 100
        self.logger.info(f"step {i}/{max_steps} ({percent:.1f}%)")
        self.logger.info("MD start")
        self.logger.info(f"doing {self.steps_per_step} MD steps")
        self.system.do_steps(self.steps_per_step)
        self.logger.info("MD finished")
        if max_steps > 0:
            ns_so_far = self.step_ns * i
            ns_total = self.step_ns * max_steps
            time_left = self.last_step_time * (max_steps - i)
            time_left_fmt = time.strftime("%H:%M:%S", time.gmtime(time_left))
            sys.stdout.write(f"\rStep {i}/{max_steps}\t"
                             f"{ns_so_far:.2f}/{ns_total:.2f}ns\t"
                             f"{percent:.1f}%\t"
                             f"ETL {time_left_fmt}\t"
                             f"{self.reactions} reactions")
        self.logger.info("D/M start")
        new_reactions = self.top.detection_modification(i)
        self.reactions += new_reactions
        self.logger.info("D/M finished")
        if new_reactions > 0:
            self.logger.info("reinitialize start")
            self.system.reinitialize()
            self.logger.info("reinitialize finished")
        if i % self.xtc_every == 0:
            self.logger.info("XTC write start")
            self.system.write_xtc_frame(self.steps_per_step)
            self.logger.info("XTC write finished")
        end_time = time.time()
        self.last_step_time = (
            0.97 * self.last_step_time + 0.03 * (end_time - start_time)
            if self.last_step_time > 0. else (end_time - start_time)
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
