#!/usr/bin/env python3
"""Very quick script to generate system.gro and expected.reactions."""

import numpy as np

from martini_daemon import PeriodicBox, write_geometry

title_line = "# frame,reaction_name;reactant1_name,reactant1_id,atoms...;reactantn_name,reactantn_id,atoms..."
reactions = []

atom_names = []
res_names = []
res_ids = []

# synchronize this with system.top manually!
reaction_name = "tri"
r1 = "x"
r2 = "y"
r3 = "z"
n_molecules = 1000
d1_range = (0, 0.3)
d2_range = (0.2, 0.4)
angle_range = (np.pi / 3, np.pi / 3 * 2)

pos = np.empty((n_molecules * 3, 3))

for i in range(n_molecules):
    # NOTE: if it starts not passing after re-generating, check whether there is an exact equality to r_max/r_min
    d1 = np.random.random() * 0.5
    d2 = np.random.random() * 0.5
    angle = np.random.random() * np.pi

    # particle 1 at x=d1, y=0
    pos[i * 3, 0] = d1
    pos[i * 3, 1] = 0

    # particle 2 at x=0, y=0
    pos[i * 3 + 1, 0] = 0
    pos[i * 3 + 1, 1] = 0

    # particle 3 on a circle with radius d2 and angle
    pos[i * 3 + 2, 0] = np.cos(angle) * d2
    pos[i * 3 + 2, 1] = np.sin(angle) * d2

    # z for all 3 same
    pos[i * 3 : i * 3 + 3, 2] = i

    # would this react?
    if (
        d1_range[0] <= d1 <= d1_range[1]
        and d2_range[0] <= d2 <= d2_range[1]
        and angle_range[0] < angle < angle_range[1]
    ):
        reactions.append(
            f"0,tri;{r1},{i * 3},{i * 3};{r2},{i * 3 + 1},{i * 3 + 1};{r3},{i * 3 + 2},{i * 3 + 2}"
        )

    atom_names.append(r1)
    atom_names.append(r2)
    atom_names.append(r3)
    res_names.append(r1)
    res_names.append(r2)
    res_names.append(r3)
    res_ids.append(i * 3 + 1)
    res_ids.append(i * 3 + 2)
    res_ids.append(i * 3 + 3)


with open("expected.reactions", "w") as f:
    f.write(f"{title_line}\n")
    f.write("\n".join(reactions))
    f.write("\n")

write_geometry(
    "system.gro",
    "example trimolecular reaction",
    atom_names,
    res_names,
    res_ids,
    PeriodicBox.cubic(n_molecules + 5.0),
    pos,
)
