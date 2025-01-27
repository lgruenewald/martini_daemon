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


class DaemonSimulation():

    system: SysStar
    top: TopStar
    gro: GromacsGroFile

    reaction_matrix: dict[(str, str), ReactionTemplate]
    initiator_list: list[Fragment]

    i: int = 0
    reactions: int = 0
    max_steps: int
    steps_per_step: int
    silent: bool  # if silent no files or stdout are written to
    traj_path: str  # trajectory to write
    out_path: str  # final geometry to write

    def __init__(self, top_path, gro_path, T=300., p=1., dt=20*femtosecond,
                 max_steps=100, steps_per_step=5000, traj_path="traj.xtc",
                 out_path="final.gro", silent=False, platform=None,
                 minimize_energy=True, generate_velocities=True,
                 remove_com_motion=True, epsilon_r=15.0,
                 nonbonded_cutoff=1.1*nanometer, include_dir=None,
                 defines={}):
        self.system, self.top = DaemonTopFile(
            top_path,
            include_dir=include_dir, defines=defines,
            epsilon_r=epsilon_r, nonbonded_cutoff=nonbonded_cutoff
        )
        self.gro = GromacsGroFile(gro_path)
        if platform is not None:
            platform = mm.Platform.getPlatformByName(platform)

        if p is not None:
            self.system.add_force(mm.MonteCarloBarostat(p, T))
        if remove_com_motion:
            self.system.add_force(mm.CMMotionRemover())
        self.reaction_matrix = self.top.build_reaction_matrix()
        self.initiator_list = self.top.get_initiator_list()

        integrator = mm.LangevinIntegrator(T, 10.0, dt)
        box = self.gro.getPeriodicBoxVectors()

        if platform is not None:
            self.system.build_context(integrator, box, platform)
        else:
            self.system.build_context(integrator, box)

        self.system.set_positions(self.gro.getPositions(True))
        if generate_velocities:
            self.system.generate_velocities(T)
        if minimize_energy:
            self.system.minimize_energy()
        if not silent:
            self.system.set_xtc_path(traj_path)
        self.max_steps = max_steps
        self.steps_per_step = steps_per_step
        self.silent = silent
        self.traj_path = traj_path
        self.out_path = out_path

    def simulate(self):
        for i in range(self.max_steps):
            self.step()
        if not self.silent:
            print()
            self.system.write_gro(self.out_path)

    def step(self):
        self.system.do_steps(self.steps_per_step)
        if not self.silent:
            sys.stdout.write(f"\rStep {self.i+1:8} of {self.max_steps}   "
                             f"[reactions: {self.reactions}]")
        self.i += 1

        pos, box = self.system.get_positions()

        pairs = []
        skip = set()
        self.top.pre_detection()
        # Detection algorithm
        for i, frag1 in enumerate(self.initiator_list):
            if frag1 is None:
                continue
            if i in skip:
                continue
            for j, frag2 in enumerate(self.initiator_list):
                if frag2 is None:
                    continue
                if j in skip:
                    continue
                rx = self.reaction_matrix.get((frag1.name, frag2.name))
                if rx is not None:
                    if self.top.detection(frag1, frag2, rx, pos, box):
                        skip.add(i)
                        skip.add(j)
                        pairs.append((frag1, frag2, rx))
                        break  # skip i - break entire loop over js with i

        if len(pairs) == 0:
            return

        # Modification algorithm
        for frag1, frag2, rx in pairs:
            self.reactions += 1
            self.top.modification(frag1, frag2, rx)

        # reinitialize context, initator list
        self.initiator_list = self.top.get_initiator_list()
        self.system.reinitialize()


if __name__ == "__main__":
    argv = sys.argv
    argc = len(sys.argv)

    if argc != 3:
        print("Usage: ./main.py <top file> <gro file>")
        quit(1)

    top_path = argv[1]
    gro_path = argv[2]

    sim = DaemonSimulation(top_path, gro_path)
    for i in range(sim.max_steps):
        sim.step()
    sys.stdout.write("\n")
    sim.system.dump()
    sim.top.dump()
    sim.system.write_gro("final.gro")
