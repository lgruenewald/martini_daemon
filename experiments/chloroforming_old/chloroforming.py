#!/usr/bin/env python3

import os
import sys
import openmm as mm
import openmm.app as mmapp
from openmm.unit import kelvin, picosecond, femtosecond,kilojoule_per_mole, \
                        kilojoule, mole, nanometer
import numpy as np
import math

sys.path.append("../../src")
from daemon_top_parser import DaemonTopFile

system, top = DaemonTopFile("system.top")
gro = mmapp.GromacsGroFile("system.gro")
system.add_force(mm.MonteCarloBarostat(1.0, 300.))
system.use_atom_type("X2")
system.build_context(mm.LangevinIntegrator(300., 1., 20 * femtosecond),
                     gro.getPeriodicBoxVectors())
system.set_positions(gro.getPositions(True))

system.set_xtc_path("traj.xtc")

for i in range(200):
    # replace 200 atoms with chloroform
    print(f"step {i}")
    for j in range(10):
        system.do_steps(500)
    # need to be updated:
    # 1. change _part_list (only affects gro writing atm)
    system._part_list[i] = ("CLF", "X2", 0., 72.)
    # 2. change charge, LJ type in _nb_force
    X2_index = system._used_atom_types["X2"]
    system._nb_force.setParticleParameters(i, [X2_index, 0.])
    # 3. add es_self_correction_forces - not done with 0 charge
#    system._es_self_correction_force.addBond(i, i, [0.5])
    # 4. change mass
    system._system.setParticleMass(i, 72.)
    # 5. reinitialize context
    system._context.reinitialize(preserveState=True)

for i in range(500):
    print(f"final eq {i}")
    # get some more frames at the end
    system.do_steps(500)

system.write_gro("final.gro")


