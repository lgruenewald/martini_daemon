#!/usr/bin/env python
# quick python script that parses Daemon .log files to generate a cost
# overview of different components (MD, D/M, reinit, xtc write)

from datetime import datetime
import sys
import os
from martini_daemon.utils import backup_try
import numpy as np

date_format = "%Y-%m-%d %H:%M:%S,%f"

clean = len(sys.argv) > 1 and sys.argv[-1] == "clean" or False
if clean:
    sys.argv.pop(-1)
path = len(sys.argv) > 1 and sys.argv[1] or "."


def extract(path):
    categories_starts = {}
    categories_sums = {}
    categories_values = {}
    md_steps = ""
    try:
        with open(path, "r") as f:
            for line in f.readlines():
                date = line[:23]
                content = line[24:]
                parsed_date = datetime.strptime(date, date_format)
                if "start" in content:
                    category = content.replace(" start", "").strip()
                    categories_starts[category] = parsed_date
                if "finished" in content:
                    category = content.replace(" finished", "").strip()
                    prev_start = categories_starts.get(category)
                    if prev_start is None:
                        continue
                    diff = parsed_date - prev_start
                    prev_diff = categories_sums.get(category) or 0.0
                    if categories_values.get(category) is None:
                        categories_values[category] = []
                    categories_sums[category] = diff.total_seconds() + prev_diff
                    categories_values[category].append(diff.total_seconds())
                if "MD steps" in content:
                    md_steps = content.replace("doing ", "").replace(" MD steps", "").strip()
    except (UnicodeDecodeError, ValueError):
        # not a text file / not a log file with dates
        return False
    if len(categories_sums) == 0:
        # not a daemon log file we recognize
        return False
    print(f"== Start of {path} ==")
    print(f"MD steps per step: {md_steps}")
    for category, seconds in categories_sums.items():
        print(f"{category}: total of {seconds} s")
        print(f"{category} n {len(categories_values[category])}")
        print(f"{category} avg {np.mean(np.array(categories_values[category]))}")
        print(f"{category} std {np.mean(np.std(categories_values[category]))}")
    print(f"== End of {path} ==")
    return True


def extract_dir(path):
    for dirpath, dirs, files in os.walk(path):
        for file in files:
            if file[0] == "#":
                continue
            filepath = os.path.join(path, dirpath, file)
            if extract(filepath) and clean:
                backup_try(filepath)


if os.path.isfile(path):
    extract(path)
else:
    extract_dir(path)
