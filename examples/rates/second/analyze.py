#!/usr/bin/env python3

from martini_daemon.helpers.rates import parse_frag_counts, parse_old_reactions

import numpy as np
import matplotlib.pyplot as plt

dm_freq = 100
xtc_freq = 5000
md_steps = 1000000

paths = []

for rate, highest_prob in [("0", 0.1), ("1", 0.1), ("1", 0.5), ("1", 1.0)]:
    for rep in range(3):
        paths.append(f"out_r{rate}_p{highest_prob}_rep{rep}")

print("Possible input files:")
for i, path in enumerate(paths):
    print(i, path)

print("choice?")
choice = input()
if len(choice) == 0:
    quit()
path = paths[int(choice)]

frag_counts = parse_frag_counts(path + ".frags", dm_freq, md_steps)
rx_counts = parse_old_reactions(path + ".reactions", dm_freq, md_steps)

x = np.linspace(dm_freq, md_steps, md_steps // dm_freq)

for name, data in frag_counts.items():
    print(f"Frag name {name} loaded")
    plt.plot(x, data[1:], label=name)

plt.legend()
plt.title("Frag counts over simulation steps")
plt.show()

for name, data in rx_counts.items():
    print(f"Reaction name {name} loaded")
    plt.plot(x, np.cumsum(data[1:]), label=name)

plt.legend()
plt.title("Rx counts over simulation steps")
plt.show()
