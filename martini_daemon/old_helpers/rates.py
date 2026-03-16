# Helpers for analyzing .frags and .reactions to produce graphs about rates

import numpy as np


def parse_old_reactions(path, interval, frame_count):
    """
    Parses the old (before June 9 2025) .reactions format,
    returns a dict of reaction name -> reaction count per D/M frame
    (np.ndarray of shape (frame_count//interval,) and type np.int64)

    args:
    path - path to the .reactions file
    interval - how often was there a D/M algo
    frame_count - how many total MD frames were written

    the latter two arguments are needed to produce the right dimension array
    even if there is no reactions in some frames.
    """

    c_frame = 0
    reactions: dict[str, np.ndarray] = {}

    with open(path, "r") as f:
        lines = f.read().split("\n")
        for line in lines:
            if line[:5] == "Frame":
                c_frame = int(line[6:].strip())
                if c_frame > frame_count:
                    raise ValueError(
                        "frame_count too low, .reactions contains frame "
                        f"{c_frame}"
                    )
            elif line[:8] == "Reaction":
                name = line.split()[1]
                if reactions.get(name) is None:
                    reactions[name] = np.zeros(
                        frame_count // interval + 1, np.int64
                    )
                reactions[name][c_frame // interval] += 1
    return reactions


def parse_frag_counts(path, interval, frame_count):
    """
    Parses the .frags format,
    returns a dict of frag name -> frag count per D/M frame
    (np.ndarray of shape (frame_count//interval,) and type np.int64)

    args:
    path - path to the .frags file
    interval - how often was there a D/M algo
    frame_count - how many total MD frames were written

    the latter two args are needed for the dimension of the array
    """

    counts: dict[str, np.ndarray] = {}

    with open(path, "r") as f:
        lines = f.read().split("\n")
        for line in lines:
            if line[:5] != "Frame":
                continue
            toks = [
                tok.strip().split(":") for tok in line.strip().split(",")
            ]
            assert toks[0][0] == "Frame"
            c_frame = int(toks[0][1])
            assert c_frame // interval <= frame_count // interval, (
                f"frame {c_frame} found in .frags, but frame_count is smaller."
                f" c_frame // interval is {c_frame // interval}"
                f" frame_count // interval is {frame_count // interval}"
            )
            for [name, count] in toks[1:]:
                if counts.get(name) is None:
                    counts[name] = np.zeros(
                        frame_count // interval + 1, np.int64
                    )
                counts[name][c_frame // interval] = int(count)
    return counts
