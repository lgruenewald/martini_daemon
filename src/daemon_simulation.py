#!/usr/bin/env python3
"""Daemon Simulation object (prototype)
the D/M algorithm + wrappers
"""

from daemon_top_parser import DaemonTopFile
from sysstar import SysStar
from topstar import TopStar, ReactionTemplate, Fragment
import sys
import openmm as mm
from openmm.app import GromacsGroFile
from openmm.unit import femtosecond
import utils


class DaemonSimulation():

    system: SysStar
    top: TopStar
    gro: GromacsGroFile

    reaction_matrix: dict[(str, str), ReactionTemplate]
    initiator_list: list[Fragment]

    def __init__(self, top_path, gro_path, T=300., p=1., dt=20*femtosecond):
        self.system, self.top = DaemonTopFile(top_path)
        self.gro = GromacsGroFile(gro_path)

        # FIXME: choice of coupling options etc
#        self.system.add_force(mm.MonteCarloBarostat(p, T))
        self.reaction_matrix, self.initiator_list = self.top.build_reaction_matrix()

        self.system.build_context(mm.LangevinIntegrator(T, 1.0, dt),
                                  self.gro.getPeriodicBoxVectors())

        self.system.set_positions(self.gro.getPositions(True))
        # FIXME: gen velocities or load velocities explicitly
        self.system.generate_velocities(T)
        # FIXME: hardcoded path
        self.system.set_xtc_path("traj.xtc")

    def step(self):
        # FIXME: hardcoded everything
        self.system.do_steps(1000)

        state = self.system.get_state()
        pos = state.getPositions()
        box = state.getPeriodicBoxVectors()[0].x

        # Detection algorithm
        for i, frag1 in enumerate(self.initiator_list):
            for j, frag2 in enumerate(self.initiator_list):
                if j >= i:
                    continue
                rx = self.reaction_matrix.get((i, j))
                if rx is not None:
                    init1 = frag1.particles[0]
                    init2 = frag2.particles[0]
                    dist = utils.pdist(pos[init1], pos[init2], box)
                    if dist < rx.max_distance:
                        print("Within cutoff distance!")

        # Modification algorithm

        # TODO destroy fragments that overlap
        # TODO destroy forces inside fragment
        # TODO instantiate product
        # TODO update initiator list
        # TODO reinitialize context

    def dump(self):
        self.system.dump()
        self.top.dump()
        self.system.write_gro("dump.gro")


if __name__ == "__main__":
    argv = sys.argv
    argc = len(sys.argv)

    if argc != 3:
        print("Usage: ./main.py <top file> <gro file>")
        quit(1)

    top_path = argv[1]
    gro_path = argv[2]

    sim = DaemonSimulation(top_path, gro_path)
    for i in range(10):
        sim.step()
    sim.dump()
