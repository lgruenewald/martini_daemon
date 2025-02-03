#!/usr/bin/env python

import os
import sys
import time
from datetime import datetime
import json
import numpy as np
import platform
import psutil
import subprocess
import openmm.version

from martini_daemon.simulation import DaemonSimulation
from martini_daemon.utils import backup_try

result_path = "benchmark_results.txt"
backup_try(result_path)
mm_platform = "CUDA"


def bprint(message, yellow=False, file_only=False):
    with open(result_path, "a") as f:
        print(f"[{datetime.now()}]", message, file=f)
    if not file_only:
        if yellow:
            print("\x1b[1;33m" + message + "\x1b[0m")
        else:
            print(message)


def test_daemon(reactive):
    defines = {}
    if reactive:
        defines["REACT"] = 1
    with open("params.json") as f:
        data = json.load(f)
    sim = DaemonSimulation(data["top"], data["gro"], max_steps=data["steps"],
                           steps_per_step=data["per_step"],
                           platform=mm_platform, defines=defines)
    sim.simulate()


def test_gromacs():
    with open("params.json") as f:
        data = json.load(f)
    top = data["top"]
    gro = data["gro"]
    steps = data["steps"] * data["per_step"]
    os.system(f"../gmxrun.sh {top} {gro} {steps}")


def bench_function(function, args, benchmark, label, n=3):
    """
        Simple helper that runs function with args n times, measures the time,
        prints the time and stddev to stdout and returns the average time.
    """

    times = []
    for i in range(n):
        start = time.time()
        function(*args)
        end = time.time()
        bprint(f"{benchmark}/{label}/{i}:{end-start:.2f}s")
        times.append(end - start)
    times = np.array(times)

    mean = np.mean(times)
    stdev = np.std(times)
    max = np.max(times)
    min = np.min(times)
    bprint(f"{benchmark}/{label} runs:{n} avg:{mean:.2f}s [{max:.2f}s..{min:.2f}s] "
           f"std:{stdev:.2f}s", yellow=True)
    return mean


def get_cpu_model():
    command = "cat /proc/cpuinfo"
    info = subprocess.check_output(command, shell=True).decode().strip()
    for line in info.split("\n"):
        if "model name" in line:
            return line.split(":")[1].strip()


def get_gpu_stats():
    command = "nvidia-smi --query-gpu=utilization.gpu,memory.total,memory.used,temperature.gpu,gpu_name --format=csv"
    info = subprocess.check_output(command, shell=True).decode().strip()
    res = ""
    for line in info.split("\n")[1:]:
        data = line.split(",")
        res += f"GPU {data[4]}: GPU utilization {data[0]}, Total mem {data[1]}, Used mem {data[2]}, Temperature {data[3]}C\n"
    res.strip()
    return res


# print data about hardware and current utilization of resources
bprint("== System information ==", file_only=True)
bprint(f"{platform.platform()}", file_only=True)  # OS info
bprint(f"CPU model: {get_cpu_model()}", file_only=True)  # CPU model
bprint(f"CPU percent used {psutil.cpu_percent()}%", file_only=True)  # current CPU utilization
mem = psutil.virtual_memory()
bprint(f"RAM total {mem.total / 1024**3:.1f}G - used {mem.percent}%", file_only=True)  # RAM amount, utilization
bprint(f"{get_gpu_stats()}", file_only=True)  # All GPU data we need
# gromacs and openmm data (TODO more?)
bprint("== Software information ==", file_only=True)
bprint(f"OpenMM version: {openmm.version.version}", file_only=True)
bprint(f"OpenMM platform being used: {mm_platform}")  # most relevant info, print it to the stdout as well

# run the benchmarks

argv = sys.argv
argc = len(argv)

benchmarks = os.listdir(os.curdir)

if argc > 1:
    tests = argv[1:]

result_path = "../" + result_path
for x in benchmarks:
    if os.path.isdir(x):
        os.chdir(x)

        bprint(f"== Benchmark {x} ==", yellow=True)
        bench_function(test_gromacs, [], x, "gromacs")
        daemon_time = bench_function(test_daemon, [True], x, "daemon with reactions")
        no_reaction_time = bench_function(test_daemon, [False], x, "daemon without reactions")

        percent = (daemon_time - no_reaction_time) / no_reaction_time * 100
        bprint(f"Reactions cost a {percent:.1f}% slowdown", yellow=True)

        os.system("../cleanup.sh")
        os.chdir("..")
