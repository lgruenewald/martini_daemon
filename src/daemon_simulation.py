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
from openmm.unit import femtosecond, nanometer
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
        self.reaction_matrix = self.top.build_reaction_matrix()
        self.initiator_list = self.top.get_initiator_list()

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
        pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)
        box = state.getPeriodicBoxVectors()[0].x

        pairs = []
        skip = set()
        # FIXME this skip system does not handle overlapping reactive fragments
        # Detection algorithm
        for i, frag1 in enumerate(self.initiator_list):
            if frag1 is None:
                continue
            for j, frag2 in enumerate(self.initiator_list):
                if frag2 is None:
                    continue
                if j >= i:
                    continue
                if i in skip or j in skip:
                    continue
                rx = self.reaction_matrix.get((frag1.name, frag2.name))
                if rx is not None:
                    if self.top.detection(frag1, frag2, rx, pos, box):
                        skip.add(i)
                        skip.add(j)
                        pairs.append((frag1, frag2, rx))

        if len(pairs) == 0:
            return

        # Modification algorithm
        for frag1, frag2, rx in pairs:
            # TODO constraint checking here
            self.top.modification(frag1, frag2, rx.p1)

        print("Modification algo ran")
        self.system.dump()
        self.top.dump()

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
    max_steps = 1000
    for i in range(max_steps):
        sim.step()
        sys.stdout.write(f"\rStep {i+1} of {max_steps}          ")
    sys.stdout.write("\n")
    sim.system.write_gro("final.gro")
