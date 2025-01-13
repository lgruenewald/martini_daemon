#!/usr/bin/env python

import os
import sys
import time

from martini_daemon.simulation import DaemonSimulation


def test_daemon(top, gro, reactive, steps, per_step):
    start = time.time()
    sim = DaemonSimulation(top, gro, max_steps=steps, steps_per_step=per_step,
                           silent=True)
    sim.simulate()
    end = time.time()
    return end - start


argv = sys.argv
argc = len(argv)

benchmarks = os.listdir(os.curdir)

if argc > 1:
    tests = argv[1:]

for x in benchmarks:
    if os.path.isdir(x):
        os.chdir(x)

        steps = 100
        per_step = 5000

        print(f"Benchmark {x}:")
        daemon_time = test_daemon(
            "system.top", "system.gro", True, steps, per_step
        )
        print(f"Reactions on: {daemon_time} seconds")
        no_reaction_time = test_daemon(
            "no_react.top", "system.gro", False, steps, per_step
        )

        print(f"Reactions off {no_reaction_time} seconds")
        percent = daemon_time / no_reaction_time * 100
        print(f"Reactions cost a {percent:.2f}% slowdown")
        print(f"({steps} reaction steps, {per_step} steps per reaction step)")

        os.chdir("..")
