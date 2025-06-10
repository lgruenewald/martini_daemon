#!/usr/bin/env python
# quick python script that parses Daemon .log files to generate a cost
# overview of different components (MD, D/M, reinit, xtc write)

from datetime import datetime
import sys

date_format = "%Y-%m-%d %H:%M:%S,%f"

path = sys.argv[1]


def extract(path):
    categories_starts = {}
    categories_sums = {}
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
                    categories_sums[category] = diff.total_seconds() + prev_diff
    except (UnicodeDecodeError, ValueError):
        # not a text file / not a log file with dates
        return False
    if len(categories_sums) == 0:
        # not a daemon log file we recognize
        return False
    print(f"== Start of {path} ==")
    for category, seconds in categories_sums.items():
        print(f"{category}: total of {seconds} s")
    print(f"== End of {path} ==")
    return True


extract(path)
