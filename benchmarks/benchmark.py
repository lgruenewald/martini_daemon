#!/usr/bin/env python

import os
import sys
import time
import json
import numpy as np

from martini_daemon.simulation import DaemonSimulation


def test_daemon(reactive):
    defines = {}
    if reactive:
        defines["REACT"] = 1
    with open("params.json") as f:
        data = json.load(f)
    sim = DaemonSimulation(data["top"], data["gro"], max_steps=data["steps"],
                           steps_per_step=data["per_step"],
                           platform="CUDA", defines=defines)
    sim.simulate()


def test_gromacs():
    with open("params.json") as f:
        data = json.load(f)
    top = data["top"]
    gro = data["gro"]
    steps = data["steps"] * data["per_step"]
    os.system(f"../gmxrun.sh {top} {gro} {steps}")


def bench_function(function, args, label, n=3):
    """
        Simple helper that runs function with args n times, measures the time,
        prints the time and stddev to stdout and returns the average time.
    """

    times = []
    for i in range(n):
        start = time.time()
        function(*args)
        end = time.time()
        print(f"{label} run {i} took {end-start:.2f}s")
        times.append(end - start)
    times = np.array(times)

    mean = np.mean(times)
    stdev = np.std(times)
    max = np.max(times)
    min = np.min(times)
    print("\x1b[1;33m"  # yellow
          f"Benchmark {label} avg:{mean:.2f}s [{max:.2f}s..{min:.2f}s] "
          f"std:{stdev:.2f}s"
          "\x1b[0m")  # reset format


argv = sys.argv
argc = len(argv)

benchmarks = os.listdir(os.curdir)

if argc > 1:
    tests = argv[1:]

for x in benchmarks:
    if os.path.isdir(x):
        os.chdir(x)

        print(f"== Benchmark {x} ==")
        bench_function(test_gromacs, [], "gromacs")
        daemon_time = bench_function(test_daemon, [True], "daemon with reactions")
        no_reaction_time = bench_function(test_daemon, [False], "daemon without reactions")

        percent = (daemon_time - no_reaction_time) / no_reaction_time * 100
        print("\x1b[1;33m"
              f"Reactions cost a {percent:.2f}% slowdown"
              "\x1b[0m")

        os.chdir("..")
