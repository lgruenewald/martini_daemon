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
from datetime import datetime
from .utils import backup_try
from math import isclose
import logging


class DaemonSimulation():

    system: SysStar
    top: TopStar
    gro: GromacsGroFile

    reaction_list: list[ReactionTemplate]
    initiator_list: dict[str, list[Fragment]]

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
        self.log_path = log_path
        backup_try(log_path)
        logging.basicConfig(filename=log_path, level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        self.log("__init__ in DaemonSimulation")
        self.log("Parsing start")
        self.system, self.top = DaemonTopFile(
            top_path,
            include_dir=include_dir, defines=defines,
            epsilon_r=epsilon_r, nonbonded_cutoff=nonbonded_cutoff
        )
        self.log("Parsing done")
        self.gro = GromacsGroFile(gro_path)
        if platform is not None:
            platform = mm.Platform.getPlatformByName(platform)

        if p is not None:
            self.system.add_force(mm.MonteCarloBarostat(p, T))
        if remove_com_motion:
            self.system.add_force(mm.CMMotionRemover())
        self.log("Building reaction matrix and initiator list")
        self.reaction_list = self.top.build_reaction_list()
        self.initiator_list = self.top.get_initiator_list()
        self.log("Reaction matrix and initiator list built")

        integrator = mm.LangevinIntegrator(T, 10.0, dt)
        box = self.gro.getPeriodicBoxVectors()

        self.log("Building context")
        if platform is not None:
            self.system.build_context(integrator, box, platform)
        else:
            self.system.build_context(integrator, box)
        self.log("Context built")

        self.log("Setting positions")
        self.system.set_positions(self.gro.getPositions(True))
        for rep in reporters:
            self.system.add_reporter(rep)
        # must set xtc path after adding reporters currently
        self.system.set_xtc_path(traj_path)
        self.system.write_xtc_frame()
        if generate_velocities:
            self.log("Generating velocities")
            self.system.generate_velocities(T)
        if minimize_energy:
            self.log("Minimizing energy")
            self.system.minimize_energy()
        self.system.write_xtc_frame()
        self.max_steps = max_steps
        self.steps_per_step = steps_per_step
        self.traj_path = traj_path
        self.out_path = out_path
        self.log("__init__ finished")

    def simulate(self):
        self.i = 0
        for i in range(self.max_steps):
            self.step()
        print()
        self.system.write_gro(self.out_path)

    def step(self):
        self.log(f"step {self.i}, doing MD steps")
        self.system.do_steps(self.steps_per_step)
        self.log(f"{self.steps_per_step} MD steps performed")
        sys.stdout.write(f"\rStep {self.i+1:8} of {self.max_steps}   "
                         f"[reactions: {self.reactions}]")
        self.i += 1

        self.log("Detection start")
        pos, box = self.system.get_positions()

        reactions = []
        skip = set()
        self.top.pre_detection()
        # Detection algorithm, generic for all reacting molecule amounts
        for rx in self.reaction_list:
            # get a product of possible reactant combinations
            n_reactants = len(rx.reactants)
            n_types_per_reactant = []  # how many frags of such reactant are in the system
            n_types = 1
            for r in rx.reactants:
                if self.initiator_list.get(r) is None:
                    # reaction isn't possible, no reactant available
                    n_types = 0
                    break
                n = len(self.initiator_list[r])
                n_types_per_reactant.append(n)
                n_types *= n
            for i in range(n_types):
                if i > 0 and i % 5000000 == 0:
                    # logging for very slow D/M algos
                    self.log(f"D algorithm ({i/n_types*100.0:.1f}%): currently doing combination {i} out of {n_types}")
                frag_ids = []
                frags = []
                frag_names = []
                remainder = i
                cont = False
                for rid in range(n_reactants):
                    frag_id = remainder % n_types_per_reactant[rid]
                    frag = self.initiator_list[rx.reactants[rid]][frag_id]
                    frag_name = frag.name
                    remainder = remainder // n_types_per_reactant[rid]
                    # continue if skip before skip cutoff is in skip
                    if (rx.skip is None or rid < rx.skip) and (frag_name, frag_id) in skip:
                        cont = True
                        break
                    # frag id's must be in order if the name is the same to
                    # prevent double counting and self reaction
                    for prev_rid, prev_frag_id in enumerate(frag_ids):
                        # always ban self reaction
                        if rx.reactants[rid] == rx.reactants[prev_rid] and \
                                frag_id == prev_frag_id:
                            cont = True
                            break
                        # don't cross the skip boundary here
                        if rx.skip is not None and prev_rid < rx.skip and rid >= rx.skip:
                            continue
                        # order them well, if it doesn't cross a skip boundary
                        if rx.reactants[rid] == rx.reactants[prev_rid] and \
                                frag_id < prev_frag_id:
                            cont = True
                            break
                    if cont:
                        break
                    frag_ids.append(frag_id)
                    frags.append(frag)
                    frag_names.append(frag_name)
                if cont:
                    continue
                # detection
                if self.top.detection(frags, rx, pos, box):
                    # if rx.skip is defined, only run the modification on the
                    # first skip atoms, and only skip the first skip atoms
                    # the other "reactants" were there only for the detection
                    if rx.skip is not None:
                        frag_ids = frag_ids[:rx.skip]
                        frags = frags[:rx.skip]
                    for j, frag_id in enumerate(frag_ids):
                        skip.add((frag_names[j], frag_id))
                    reactions.append((frags, rx))

        self.log("Detection finished")
        if len(reactions) == 0:
            return

        self.log("Modification start")
        self.top.pre_modification()
        # Modification algorithm
        for (frags, rx) in reactions:
            self.reactions += 1
            self.top.modification(frags, rx)
        self.top.post_modification()

        # reinitialize context, initator list
        self.initiator_list = self.top.get_initiator_list()
        self.log("Modification finished")
        self.log("reinitializing")
        self.system.reinitialize()
        self.log("reinitialized")

    def log(self, message):
        self.logger.info(f"{datetime.now()} {message}")


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
