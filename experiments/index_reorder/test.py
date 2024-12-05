#!/usr/bin/env python

import openmm as mm

sys = mm.System()

f1 = mm.HarmonicBondForce()
f2 = mm.HarmonicAngleForce()

i1 = sys.addForce(f1)
i2 = sys.addForce(f2)

print(sys.getForces())
sys.removeForce(i1)

print(sys.getForces())
print(f2)
print(f2 == sys.getForce(0))
