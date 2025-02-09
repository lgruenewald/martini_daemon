#!/usr/bin/env python3
"""Daemon Simulation object (prototype)
the D/M algorithm + wrappers
"""

from .top_parser import DaemonTopFile
from .sysstar import SysStar
from .topstar import TopStar, ReactionTemplate, Fragment
import sys
import openmm as mm
from openmm.app import GromacsGroFile
from openmm.unit import femtosecond, nanometer
from .utils import backup_try
import logging


class DaemonSimulation():

    system: SysStar
    top: TopStar
    gro: GromacsGroFile

    reaction_list: list[ReactionTemplate]
    init_list: list[Fragment]
    init_map: dict[str, list[int]]

    i: int = 0
    reactions: int = 0
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
                 defines={}, log_path=None, reporters=[]):
        if traj_path is None:
            traj_path = sim_name + ".xtc"
        if out_path is None:
            out_path = sim_name + ".gro"
        if log_path is None:
            log_path = sim_name + ".log"
        backup_try(log_path)
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s %(message)s")
        fh = logging.FileHandler(log_path)
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)
        self.logger.info("__init__ in DaemonSimulation")
        self.logger.info("Parsing start")
        self.system, self.top = DaemonTopFile(
            top_path,
            include_dir=include_dir, defines=defines,
            epsilon_r=epsilon_r, nonbonded_cutoff=nonbonded_cutoff
        )
        self.logger.info("Parsing done")
        self.gro = GromacsGroFile(gro_path)
        if platform is not None:
            platform = mm.Platform.getPlatformByName(platform)

        if p is not None:
            self.system.add_force(mm.MonteCarloBarostat(p, T))
        if remove_com_motion:
            self.system.add_force(mm.CMMotionRemover())

        integrator = mm.LangevinIntegrator(T, 10.0, dt)
        box = self.gro.getPeriodicBoxVectors()

        self.logger.info("Building context")
        if platform is not None:
            self.system.build_context(integrator, box, platform)
        else:
            self.system.build_context(integrator, box)
        self.logger.info("Context built")

        self.logger.info("Setting positions")
        self.system.set_positions(self.gro.getPositions(True))
        for rep in reporters:
            rep = rep(self.system, self.top)
            self.system.add_reporter(rep)
            self.top.add_reporter(rep)
        # must set xtc path after adding reporters currently
        self.system.set_xtc_path(traj_path)
        self.system.write_xtc_frame()
        if generate_velocities:
            self.logger.info("Generating velocities")
            self.system.generate_velocities(T)
        if minimize_energy:
            self.logger.info("Minimizing energy")
            self.system.minimize_energy()
        self.system.write_xtc_frame()
        self.max_steps = max_steps
        self.steps_per_step = steps_per_step
        self.traj_path = traj_path
        self.out_path = out_path
        self.logger.info("__init__ finished")

    def simulate(self):
        self.i = 0
        for i in range(self.max_steps):
            self.step()
        print()
        self.system.write_gro(self.out_path)

    def step(self):
        self.logger.info(f"step {self.i}:")
        self.logger.info("MD start")
        self.system.do_steps(self.steps_per_step)
        self.logger.info("MD finished")
        sys.stdout.write(f"\rStep {self.i+1:8} of {self.max_steps}   "
                         f"[reactions: {self.reactions}]")
        self.i += 1
        self.logger.info("D/M start")
        self.reactions += self.top.detection_modification()
        self.logger.info("D/M finished")
        self.logger.info("reinitialize start")
        self.system.reinitialize()
        self.logger.info("reinitialize finished")


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
