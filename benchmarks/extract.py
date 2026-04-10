#!/usr/bin/env python
# quick python script that parses Daemon .log files to generate a cost
# overview of different old_components (MD, D/M, reinit, xtc write)

import sys

from martini_daemon import extract_timings_from_log

path = sys.argv[1]

categories_sums = extract_timings_from_log(path)

print(f"== Start of {path} ==")
for category, seconds in categories_sums.items():
    print(f"{category}: total of {seconds} s")
print(f"== End of {path} ==")
