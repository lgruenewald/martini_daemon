#!/usr/bin/env python3
# New benchmark, designed to be simple to run and fast to finish in contrast
# with ./benchmark.py. Only runs 25 md steps, D/M step and xtc step once.

from martini_daemon import Simulation, extract_timings_from_log, XTCReporter, ToptrajReporter
import os
import sys
import time
import platform
import datetime

print("Hostname:", platform.node())
print("Date: ", datetime.datetime.now())

tmpdir = ".temp"

if len(sys.argv) > 1:
    inputs = sys.argv[1:]
else:
    inputs = [x for x in os.listdir("inputs/") if x[-4:] == ".top"]

for inp in inputs:
    if not os.path.isdir(tmpdir):
        os.mkdir(tmpdir)
    os.chdir(tmpdir)
    print("Input:", inp)
    start = time.time()
    sim = Simulation(
        "../inputs/" + inp, "../inputs/" + inp.replace(".top", ".gro"),
        10000,
        [XTCReporter(), ToptrajReporter()],
        25, 25,
        defines={"REACT": "1"}
    )
    start = time.time()
    for i in range(10):
        sim.step(250, traj=True, dm=True)
    break
    categories = extract_timings_from_log("out.log", ignore_first=True)
    if not categories:
        print("log for ", inp, " can't be extracted!")
        continue
    for cat, ctime in categories.items():
        print(f"{cat}: {ctime} s")
    print("")

    os.chdir("..")

