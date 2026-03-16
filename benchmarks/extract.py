#!/usr/bin/env python
# quick python script that parses Daemon .log files to generate a cost
# overview of different old_components (MD, D/M, reinit, xtc write)

from martini_daemon.old_helpers.log_extract import extract
import sys

path = sys.argv[1]

categories_sums = extract(path)

print(f"== Start of {path} ==")
for category, seconds in categories_sums.items():
    print(f"{category}: total of {seconds} s")
print(f"== End of {path} ==")
