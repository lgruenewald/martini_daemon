#!/usr/bin/env python3

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
import shutil
from datetime import datetime

from martini_daemon.old_simulation import Simulation
from martini_daemon.utils import backup_try

# setup logging
result_path = f"{datetime.now()}.log"
backup_try(result_path)
mm_platform = "CUDA"
tmpdir = ".tmp"

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

# gromacs helper
cached = {}


def bench_gromacs(top, gro, steps):
    """
    In a temporary folder that it creates and cleans up after,
    it runs a gmx simulation of top, gro for steps steps,
    and returns the time it took as a tuple: (grompp time, mdrun time)

    since it'sjust a simple reference value, we run it with n=1
    because from experience, without reactions there is low variability

    we also cache the result for the same combo of parameters,
    because it's just a simple reference value
    """
    grompp_command = f"../gmx_grompp.sh {top} {gro} {steps}"
    mdrun_command = "../gmx_mdrun.sh"
    if cached.get(grompp_command) is not None:
        grompp, mdrun = cached[grompp_command]
        logger.info(f"Cached gromacs time grompp: {grompp:.2f} s; mdrun: {mdrun:.2f} s available for {grompp_command}")
        return grompp, mdrun

    os.mkdir(tmpdir)
    os.chdir(tmpdir)
    logger.info(grompp_command)
    start_grompp = time.time()
    assert os.system(grompp_command) == 0
    end_grompp = time.time()
    grompp = end_grompp - start_grompp
    start_mdrun = time.time()
    logger.info(mdrun_command)
    assert os.system(mdrun_command) == 0
    end_mdrun = time.time()
    mdrun = end_mdrun - start_mdrun
    logger.info(f"Gromacs finished in: grompp: {grompp:.2f} s; {mdrun:.2f} s.")
    cached[grompp_command] = (grompp, mdrun)
    os.chdir("..")
    shutil.rmtree(tmpdir)
    return grompp, mdrun


# daemon n repetition setup
def setup_daemon():
    """
    called before n repetitions

    side effects:
    makes and enters a temporary folder
    """
    os.mkdir(tmpdir)
    os.chdir(tmpdir)


# daemon run
def bench_daemon(top, gro, steps, freq, reactive, force_reinit):
    """
    single repetition of daemon, returns the simulation log file name

    side effect: creation of the log file the name of which is returned
    """
    defines = {}
    if reactive:
        defines["REACT"] = "1.0"

    sim_name = f"out_{datetime.now()}"
    logger.info(
        f"Running daemon with {top} {gro} for {steps} steps {freq} "
        f"freq {reactive} reactive force_reinit {force_reinit}."
    )
    start_grompp = time.time()
    sim = Simulation(
        top_path=top, gro_path=gro, md_steps=steps,
        dm_frequency=freq, xtc_frequency=5000,
        sim_name=sim_name,
        platform=mm_platform, defines=defines,
        force_reinitialize=force_reinit
    )
    sim.minimize_energy()
    sim.generate_velocities(300)
    end_grompp = time.time()
    grompp = end_grompp - start_grompp

    logger.info(f"Parsing done in {grompp:.1f} s.")
    start_mdrun = time.time()
    sim.simulate()
    end_mdrun = time.time()
    mdrun = end_mdrun - start_mdrun
    logger.info(f"Simulation done in {mdrun:.1f} s.")

    probe_start = time.time()
    sim.step(100000, xtc=True, dm=False, neighbor=False)
    probe_end = time.time()
    probe = probe_end - probe_start
    logger.info(f"Probe 100k steps took {probe:.2f} s.")

    for h in sim.logger.handlers:
        h.flush()
        h.close()

    return (sim_name + ".log", grompp, mdrun, probe)


# log extractor, dat writer
# date format for parsing .log files
date_format = "%Y-%m-%d %H:%M:%S,%f"
# hardcoded categories we care about
categories = [
    "Parsing", "MD", "Detection", "Modification", "Reinitialize", "XTC write"
]


def extract(log_path, result_prefix, result_hist_prefix):
    """
    Takes daemon .log file at log_path, returns data about how much time
    the different old_components took (md, detection, reinit, ...)

    Also writes the processed log file data to result_paths
    and histograms to result_hist_path, in .dat format
    (example dat filename: result_prefix + <category> + ".dat")

    In the comments of each dat file the main log timestamp is given

    Returns the information required for the summary
    """
    logger.info(
        f"Extracting log info from {log_path} to "
        f"{result_prefix}* and {result_hist_prefix}*."
    )
    categories_starts = {}
    categories_values = {}
    reactions = []
    for cat in categories:
        categories_values[cat] = [0.0]
    i = 0

    with open(log_path, "r") as f:
        for line in f.readlines():
            date = line[:23]
            content = line[24:]
            parsed_date = datetime.strptime(date, date_format)
            if "step " in content:
                for cat in categories:
                    categories_values[cat].append(0.0)
                i += 1
                categories_starts = {}

            if "start" in content:
                category = content.replace(" start", "").strip()
                if category not in categories:
                    continue
                categories_starts[category] = parsed_date
            if "finished" in content:
                category = content.replace(" finished", "").strip()
                if category not in categories:
                    continue
                prev_start = categories_starts[category]
                diff = parsed_date - prev_start
                categories_values[category][i] = diff.total_seconds()
            if "reactions" in content:
                count = int(content.replace("reactions", "").strip())
                reactions.append(count)
    if len(categories_values) == 0:
        # not a daemon log file we recognize
        raise ValueError("Attempt to extract a non daemon log file")
    reactions_path = result_prefix + "reactions.dat"
    reac_arr = np.array(reactions)
    np.savetxt(
        reactions_path,
        reac_arr,
        header="#Autogenerated by benchmark.py\n"
               f"# Associated benchmark log: {result_path}\n"
               f"@reaction count",
        comments="",
        delimiter=","
    )
    sums = {}
    for cat, x in categories_values.items():
        arr = np.array(x)
        hist, edges = np.histogram(arr, bins=50)
        assert len(hist) + 1 == len(edges)
        dat_path = result_prefix + cat + ".dat"
        hist_path = result_hist_prefix + cat + ".dat"
        np.savetxt(
            dat_path,
            arr,
            header="# Autogenerated by benchmark.py\n"
            f"# Associated benchmark log: {result_path}\n"
            f"@Time (s)\n",
            comments="",
            delimiter=","
        )
        midpoints = np.zeros((len(hist),))
        for i in range(len(hist)):
            midpoints[i] = 0.5 * edges[i] + 0.5 * edges[i+1]
        np.savetxt(
            hist_path,
            np.column_stack((midpoints, hist)),
            header="# Autogenerated by benchmark.py\n"
            f"# Associated benchmark log: {result_path}\n"
            f"@midpoint,count\n",
            comments="",
            delimiter=","
        )
        sums[cat] = np.sum(arr)
    return sums


# n repetition collector
def finish_daemon(result_dir, name, reactive_runs, non_reactive_runs):
    """
    called after n repetitions, with the result of every run of bench_daemon
    as a list

    args: 
    name - name of the run, used to generate the result .dat path
    two lists - reactive and non reactive log file names

    side effects:
    generates the .dat files in results/
    <sim name>_n_*.dat and <sim name>_n_hist_*.dat:
    - The first gives the time * (md, detection, ...) took, n represents which
        repetition it is
    - The second gives the histogram for the first piece of data

    <sim name>_reactions_n.dat and <sim name>_reactions_n_hist.dat file:
    Similar to previous, but it counts the reactions, rather than telling the
    time for each component

    exits and cleans up the temporary folder, deleting the logs, xtc, ...

    returns:
    a dictionary that contains the averaged-over-repetitions
    times that each component took during a single run
    (required for the summary)
    """
    logger.info(f"Finishing runs for {name}.")
    os.chdir("..")

    reactive_data = []
    non_reactive_data = []
    reactive_parsetime = 0.
    reactive_runtime = 0.
    non_reactive_parsetime = 0.
    non_reactive_runtime = 0.
    reactive_probe = 0.
    non_reactive_probe = 0.
    for i, (run, parsetime, runtime, probe) in enumerate(reactive_runs):
        log_path = os.path.join(tmpdir, run)
        dat_name = f"{name}_{i}_"
        dat_path = os.path.join(result_dir, dat_name)
        hist_name = f"{name}_{i}_hist_"
        hist_path = os.path.join(result_dir, hist_name)

        reactive_data.append(extract(log_path, dat_path, hist_path))
        # averages over n
        reactive_parsetime = parsetime if i == 0 else (
            reactive_parsetime * (i-1) + parsetime
        ) / i
        reactive_runtime = runtime if i == 0 else (
            reactive_runtime * (i-1) + runtime
        ) / i
        reactive_probe = probe if i == 0 else (
            reactive_probe * (i-1) + probe
        ) / i

    for i, (run, parsetime, runtime, probe) in enumerate(non_reactive_runs):
        log_path = os.path.join(tmpdir, run)
        dat_name = f"{name}_NO_REACTION_{i}_"
        dat_path = os.path.join(result_dir, dat_name)
        hist_name = f"{name}_NO_REACTION_{i}_hist_"
        hist_path = os.path.join(result_dir, hist_name)

        non_reactive_data.append(extract(log_path, dat_path, hist_path))
        # averages over n
        non_reactive_parsetime = parsetime if i == 0 else (
            non_reactive_parsetime * (i-1) + parsetime
        ) / i
        non_reactive_runtime = runtime if i == 0 else (
            non_reactive_runtime * (i-1) + runtime
        ) / i 
        non_reactive_probe = probe if i == 0 else (
            non_reactive_probe * (i-1) + probe
        ) / i

    parse_slowdown = (
        (reactive_parsetime - non_reactive_parsetime)
        / non_reactive_parsetime * 100.
    )
    run_slowdown = (
        (reactive_runtime - non_reactive_runtime) / non_reactive_runtime * 100.
    )
    logger.info(
        f"Reactive systems for {name} took on average "
        f"{reactive_parsetime:.1f} s to parse {reactive_runtime:.1f} to run."
    )
    logger.info(
        f"Not reactive systems for {name} took on average "
        f"{non_reactive_parsetime:.1f} s to parse "
        f"{non_reactive_runtime:.1f} to run."
    )
    logger.info(
        f"Parse slowdown of {parse_slowdown:.1f}%. "
        f"Run slowdown of {run_slowdown:.1f}%."
    )
    logger.info(
        f"Reactive probe {reactive_probe:.2f}. "
        f"Non reactive probe {non_reactive_probe:.2f}."
    )

    # TODO also copy over the energy data (T, box size, etc.)
    shutil.rmtree(tmpdir)
    return (
        reactive_data, non_reactive_data, parse_slowdown,
        run_slowdown, reactive_probe, non_reactive_probe
    )


# main helper
def run_recipe(recipe_name, data):
    """
    args: recipe name, data json object, which is a list of dictionaries.

    A recipe is a list of reactive MD simulation parameters.

    To run the recipe, every simulation in the list
    is ran in gromacs (without reactions),
    daemon (without reactions) and daemon (with reactions).

    Afterwards, performance data is extracted to results/recipe_name.
    Every run will have various .dat files generated there.

    Additionally, generates the following summary .dat file
    in results containing the main information we care about.

    Folder structure:
    results/<recipe name>/*.dat

    Comment info:
    - timestamp
    - % slowdown for reactions, pessimistic formula

    Content of summary.dat file:
    - reactive, non reactive and gromacs for every entry
    - Plottable graph of time of each component, in seconds
    """

    logger.info(f"Starting recipe {recipe_name}.")
    result_dir = os.path.join("results", recipe_name)
    os.mkdir(result_dir)
    entries = []

    for entry in data:
        name = entry["name"]
        logger.info(f"Starting system {name}.")
        # top, gro are relative to tmpdir
        top = os.path.join("../inputs", entry["top"])
        gro = os.path.join("../inputs", entry["gro"])
        steps = entry["md_steps"]
        freq = entry["freq"]
        n = entry["n"]
        force_reinit = entry.get("force_reinitialize")

        grompp, mdrun = bench_gromacs(top, gro, steps)  # enters then exits tmp dir

        reactive_runs = []
        non_reactive_runs = []
        setup_daemon()  # enters tmp dir
        for i in range(n):
            reactive_runs.append(
                bench_daemon(top, gro, steps, freq, True, force_reinit)
            )
            non_reactive_runs.append(
                bench_daemon(top, gro, steps, freq, False, force_reinit)
            )

        # this exits the tmpdir
        (
            reac_data, non_reac_data, parse_slowdown, run_slowdown,
            rprobe, nprobe
        ) = finish_daemon(result_dir, name, reactive_runs, non_reactive_runs)
        entries.append((
            entry, reac_data, non_reac_data, grompp, mdrun, parse_slowdown,
            run_slowdown, rprobe, nprobe
        ))

    # generate summary
    summary = os.path.join(result_dir, "summary.dat")
    with open(summary, "w") as f:

        cats = ",".join(categories)
        f.write(
            "# Autogenerated by benchmark.py\n"
            f"# Associated benchmark log: {result_path}\n"
            f"@ {cats}\n"
        )

        for (
            entry, reac_data, non_reac_data, grompp, mdrun, parse_slowdown,
            run_slowdown, rprobe, nprobe
        ) in entries:
            top = entry["top"]
            gro = entry["gro"]
            steps = entry["md_steps"]
            freq = entry["freq"]
            n = entry["n"]
            name = entry["name"]
            no_reac_name = name + " no reactions"
            gro_name = name + " gromacs"
            values = []
            no_reac_values = []
            gromacs_values = []
            for cat in categories:
                # react, average of n
                sum = 0
                for x in reac_data:
                    sum += x[cat] / len(reac_data)
                values.append(f"{sum}")
                # no react, average of n
                sum = 0
                for x in non_reac_data:
                    sum += x[cat] / len(non_reac_data)
                no_reac_values.append(f"{sum}")
                # gromacs, n == 1
                if cat == "Parsing":
                    gromacs_values.append(f"{grompp}")
                elif cat == "MD":
                    gromacs_values.append(f"{mdrun}")
                else:
                    gromacs_values.append("0")

            values_str = ",".join(values)
            no_reac_str = ",".join(no_reac_values)
            gromacs_str = ",".join(gromacs_values)
            f.write(
                f"# name: {name} top: {top} gro: {gro} steps: {steps} "
                f"freq: {freq} n: {n}\n"
                f"# parse slowdown: {parse_slowdown:.1f}% "
                f"run slowdown: {run_slowdown:.1f}%\n"
                f"# reactive probe {rprobe:.2f} s; "
                f"non reactive probe {nprobe:.2f} s.\n"
                f"{name},{values_str}\n"
                f"{no_reac_name},{no_reac_str}\n"
                f"{gro_name},{gromacs_str}\n"
            )


# print data about hardware and current utilization of resources
def get_cpu_model():
    command = "cat /proc/cpuinfo"
    info = subprocess.check_output(command, shell=True).decode().strip()
    for line in info.split("\n"):
        if "model name" in line:
            return line.split(":")[1].strip()


def get_gpu_stats():
    command = (
        "nvidia-smi --query-gpu=utilization.gpu,memory.total,"
        "memory.used,temperature.gpu,gpu_name --format=csv"
    )
    info = subprocess.check_output(command, shell=True).decode().strip()
    res = ""
    for line in info.split("\n")[1:]:
        data = line.split(",")
        res += (
            f"GPU {data[4]}: GPU utilization {data[0]}, "
            f"Total mem {data[1]}, Used mem {data[2]}, "
            f"Temperature {data[3]}C\n"
        )
    res.strip()
    return res


logger.debug(f"Log file: {result_path}")
logger.debug(f"Hostname: {platform.node()}")
logger.debug(f"Platform: {platform.platform()}")
logger.debug(
    f"CPU model: {get_cpu_model()} utilization {psutil.cpu_percent()}%"
)
logger.debug(f"{get_gpu_stats()}")
logger.debug(f"OpenMM version: {openmm.version.version}")
logger.debug(f"OpenMM platform being used: {mm_platform}")
logger.debug("")

# main code: collect and run recipes
if not os.path.exists(os.path.join(os.curdir, "benchmark.py")):
    print("Error: must run benchmarks from the folder where benchmark.py is.")
    exit(1)

results_path = os.path.join(os.curdir, "results")
if not os.path.exists(results_path):
    os.mkdir(results_path)

for x in os.listdir(results_path):
    print(
        "Results folder (./result) is not empty. "
        "Please back up previous benchmark results, "
        "delete existing results folder? (y/N)"
    )
    answer = input()
    if answer in {"y", "Y"}:
        shutil.rmtree(results_path)
        os.mkdir(results_path)
        break
    else:
        print("exiting")
        exit(1)

if os.path.exists(tmpdir):
    shutil.rmtree(tmpdir)

recipes_path = os.path.join(os.curdir, "recipes")
recipes = os.listdir(recipes_path)

for recipe in recipes:
    logger.info(f"Doing recipe {recipe}.")

    with open(os.path.join(recipes_path, recipe)) as f:
        run_recipe(os.path.splitext(recipe)[0], json.load(f))

fh.flush()
fh.close()
shutil.copy(result_path, "results/benchmark_info.txt")
