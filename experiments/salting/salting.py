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
system.use_atom_type("TQ5")
system.build_context(mm.LangevinIntegrator(300., 10., 20 * femtosecond),
                     gro.getPeriodicBoxVectors())
system.set_positions(gro.getPositions(True))

system.set_xtc_path("traj.xtc")

for i in range(400):
    print(f"step {i}")
    for j in range(10):
        system.do_steps(500)
    # every time replace 2 waters with 1 na and 1 cl
    # need to be updated:
    # 1. change _part_list (only affects gro writing atm)
    is_sodium = i % 2 == 0
    part_name = "NA" if is_sodium else "CL"
    charge = 1.0 if is_sodium else -1.0
    system._part_list[i] = (part_name, "TQ5", charge, 36.)
    # 2. change charge, LJ type in _nb_force
    TQ5_index = system._used_atom_types["TQ5"]
    system._nb_force.setParticleParameters(i, [TQ5_index, charge])
    # 3. add es_self_correction_forces
    system._es_self_correction_force.addBond(i, i, [0.5])
    # 4. change mass to 36
    system._system.setParticleMass(i, 36.)
    # 5. reinitialize context
    system._context.reinitialize(preserveState=True)

for i in range(500):
    print(f"final eq {i}")
    # get 50 more frames with only salt
    system.do_steps(500)
    system._context.reinitialize(preserveState=True)

system.write_gro("final.gro")


