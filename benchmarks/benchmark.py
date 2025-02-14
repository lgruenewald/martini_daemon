#!/usr/bin/env python

import os
import sys
import time
import json
import numpy as np
import platform
import psutil
import subprocess
import openmm.version
import logging
from datetime import datetime

from martini_daemon.simulation import DaemonSimulation
from martini_daemon.utils import backup_try

result_path = "benchmark_results.log"
backup_try(result_path)
mm_platform = "CUDA"


logger = logging.getLogger(__name__)
time_formatter = logging.Formatter("%(asctime)s %(message)s")
message_only = logging.Formatter("%(message)s")
fh = logging.FileHandler(result_path)
fh.setLevel(logging.DEBUG)
fh.setFormatter(time_formatter)
sh = logging.StreamHandler(sys.stdout)
sh.setLevel(logging.INFO)
sh.setFormatter(message_only)
logger.setLevel(logging.DEBUG)
logger.addHandler(fh)
logger.addHandler(sh)


def test_daemon(reactive, data):
    defines = {}
    if reactive:
        defines["REACT"] = 1

    p = data.get("p") or 1.0
    if p == 0 or p == 0.:
        p = None
    sim = DaemonSimulation(data["top"], data["gro"], max_steps=data["steps"],
                           steps_per_step=data["per_step"], p=p,
                           platform=mm_platform, defines=defines)
    sim.simulate()
    for h in sim.logger.handlers:
        h.flush()
        h.close()
    os.rename("out.log", f"out_{datetime.now()}.log".replace(" ", "_"))


def test_gromacs(data):
    top = data["top"]
    gro = data["gro"]
    steps = data["steps"] * data["per_step"]
    p = data.get("p") or 1.0
    ptype = "C-rescale" if p > 0 else "no"
    os.system(f"../gmxrun.sh {top} {gro} {steps} {ptype} {p}")


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
        logger.info(f"{benchmark} {label} run {i}: {end-start:.2f} s")
        times.append(end - start)
    times = np.array(times)

    mean = np.mean(times)
    stddev = np.std(times)
    max = np.max(times)
    min = np.min(times)
    logger.info(f"{benchmark} {label} runs: {n} avg: {mean:.2f} s "
                f"[{max:.2f} s .. {min:.2f} s] "
                f"std dev: {stddev:.2f} s")
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
logger.debug("== System information ==")
logger.debug(f"{platform.platform()}")
logger.debug(f"CPU model: {get_cpu_model()} utilization {psutil.cpu_percent()}%")
logger.debug(f"{get_gpu_stats()}")
logger.debug("== Software information ==")
logger.debug(f"OpenMM version: {openmm.version.version}")
logger.debug(f"OpenMM platform being used: {mm_platform}")
logger.debug("")

# run the benchmarks

argv = sys.argv
argc = len(argv)

benchmarks = os.listdir(os.curdir)

if argc > 1:
    benchmarks = argv[1:]

for x in benchmarks:
    if os.path.isdir(x):
        os.chdir(x)

        with open("params.json") as f:
            data = json.load(f)

        for bench in data["benchmarks"]:
            name = bench["name"]
            n = bench["n"]

            logger.info(f"Benchmark {x}")
            logger.info(f"Parameters: {bench}")
            daemon_time = bench_function(test_daemon, [True, bench], name, "daemon with reactions", n=n)
            no_reaction_time = bench_function(test_daemon, [False, bench], name, "daemon without reactions", n=n)
            bench_function(test_gromacs, [bench], name, "gromacs", n=n)

            percent = (daemon_time - no_reaction_time) / no_reaction_time * 100
            logger.info(f"Reactions cost a {percent:.1f}% slowdown")

        os.system("../cleanup.sh")
        os.chdir("..")
