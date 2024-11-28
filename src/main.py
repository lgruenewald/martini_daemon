#!/usr/bin/env python3

from daemon_top_parser import DaemonTopFile
import utils
import sys
import openmm as mm
import openmm.app as mmapp
from openmm.unit import kelvin, picosecond, femtosecond


argv = sys.argv
argc = len(sys.argv)

if argc != 3:
    print("Usage: ./main.py <top file> <gro file>")
    quit(1)

top_path = argv[1]
gro_path = argv[2]

system, top = DaemonTopFile(top_path)


gro = mmapp.GromacsGroFile(gro_path)

mmtopol = mmapp.Topology()
mmtopol._numAtoms = system.len_particles()
mmtopol._periodicBoxVectors = gro.getPeriodicBoxVectors()

system.add_force(mm.MonteCarloBarostat(1.0, 300 * kelvin))

system.build_context(mm.LangevinIntegrator(300 * kelvin, 10.0 / picosecond,
                                           20 * femtosecond),
                     gro.getPeriodicBoxVectors())

pos = system.set_positions(gro.getPositions(True))
system.set_xtc_path("traj.xtc")
system.minimize_energy(10, 0)
system.generate_velocities(300 * kelvin)

for i in range(1000):
    system.do_steps(5)

system.write_gro("final.gro")
