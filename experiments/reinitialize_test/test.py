#!/usr/bin/env python

import openmm as mm
import openmm.app as mmapp

from openmm import Vec3
import time
import numpy as np
import os

sys = mm.System()
for i in range(1000):
    sys.addParticle(10)


f1 = mm.HarmonicBondForce()
f2 = mm.HarmonicBondForce()
f3 = mm.HarmonicBondForce()

for i, j in zip(range(100), range(100, 200)):
    f1.addBond(i, j, 1.0, 15000)


for i, j in zip(range(300, 400), range(400, 500)):
    f2.addBond(i, j, 1.0, 15000)


for i, j in zip(range(600, 700), range(700, 800)):
    f3.addBond(i, j, 1.0, 15000)

sys.addForce(f1)
sys.addForce(f2)
sys.addForce(f3)

positions = np.zeros((1000, 3))

for i in range(1000):
    positions[i] = np.array([i % 100., i / 100., 0.0])

box = [Vec3(120., 0., 0.), Vec3(0., 20., 0.), Vec3(0., 0., 20.)]


sys.setDefaultPeriodicBoxVectors(*box)

context = mm.Context(sys, mm.LangevinIntegrator(300., 10., 0.02))
context.setPositions(positions)

top = mmapp.Topology()
top._numAtoms = 1000
os.remove("out.xtc")
xtc = mmapp.XTCFile("out.xtc", top, 0.02)


def write_xtc():
    pos = context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(mm.unit.nanometer)
    xtc.writeModel(pos)


write_xtc()


def timeit(func, n=100, prep=None):
    res = np.zeros(n)
    for i in range(n):
        if prep is not None:
            prep()
        start = time.time()
        func()
        end = time.time()
        res[i] = end - start
    return np.mean(res), np.std(res)

mm.LocalEnergyMinimizer.minimize(context)
print("context set up, energy minimized")

print("timing 100 steps only")


def step100():
    context.getIntegrator().step(100)
    write_xtc()


print("step100", timeit(step100, n=100))

print("timing empty reinitialize")


def reinit():
    context.reinitialize(preserveState=True)
    context.getIntegrator().step(100)
    write_xtc()


print("step100reinit", timeit(reinit, n=100))


print("timing addBond to f1 + 100 steps")

bond = 0


def addbond():
    global bond
    f1.addBond(100+bond, 200+bond, 1.0, 100.)
    bond += 1


print("addbondreinit", timeit(reinit, prep=addbond, n=100))


print("remove bond reinit from f3 + 100 steps")

i = 100


def removebond():
    global i, f3
    sys.removeForce(2)
    f3 = mm.HarmonicBondForce()
    for j, k in zip(range(600, 600 + i), range(700, 700+i)):
        f3.addBond(j, k, 1.0, 15000)
    i -= 1
    sys.addForce(f3)


print("removebondreinit", timeit(reinit, prep=removebond, n=100))
