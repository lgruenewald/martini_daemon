"""
Readers for old formats used during Martini Daemon development, which have been replaced
by a newer format.
"""
import numpy as np
import zlib
import molly
import struct
import re

def read_compressed(path):
    with open(path, "rb") as f:
        return zlib.decompress(f.read())

def read_atoms(path):
    """
    file format - zlib compressed following:
    header: n_atoms (64 bit unsigned int)
    body - for each frame, without any padding or separator:
        - n_atoms s8 atom names
        - n_atoms s8 atom types
        - n_atoms f32 charges
        - n_atoms f32 masses
    f32 being a 32 bit float, s8 being an utf-8 string cut off at 8 bytes,
    null terminated if shorter than 8 bytes.

    note: endianness is the native endianness

    Returns:
    n_frames, n_atoms, names, types, charges, masses
    where the last 4 are numpy arrays of dimension (n_frames, n_atoms)
    """
    data = read_compressed(path)
    n_atoms, = struct.unpack("=Q", data[0:8])
    # length for a single atom: 8+8+4+4 = 24 bytes
    n_frames = (len(data) - 8) // (n_atoms * 24)
    assert (len(data) - 8) % (n_atoms * 24) == 0, "The .atoms file is incomplete or the wrong file type."
    names = np.empty((n_frames, n_atoms), dtype="a8")
    types = np.empty((n_frames, n_atoms), dtype="a8")
    charges = np.empty((n_frames, n_atoms), dtype=np.float32)
    masses = np.empty((n_frames, n_atoms), dtype=np.float32)
    for frame in range(n_frames):
        start = 8 + (frame * 24 * n_atoms)
        names[frame] = np.frombuffer(
            data[start:start+(8*n_atoms)], dtype="a8"
        )
        types[frame] = np.frombuffer(
            data[start+8*n_atoms:start+16*n_atoms], dtype="a8"
        )
        charges[frame] = np.frombuffer(
            data[start+16*n_atoms:start+20*n_atoms], dtype=np.float32
        )
        masses[frame] = np.frombuffer(
            data[start+20*n_atoms:start+24*n_atoms], dtype=np.float32
        )
    return n_frames, n_atoms, names, types, charges, masses



def read_bonds(path: str, read_n_atoms=True) -> tuple[int, int, list[np.ndarray]]:
    """
    0.1 .bonds file reader.

    File format - zlib compressed binary data:
    - n_atoms - 8 byte integer
    - bond info - for each frame:
        - n_bonds - number of bonds in this frame (64 bit unsigned integer)
        - i, j - atom indices for a bond (both 32 bit unsigned integers)
    Note: assumes at most 2^32 atoms. Endianness used is native.

    Returns: n_frames, n_atoms, frames
    frames = a list of frames, each frame containing a numpy
    array of (n_bonds, 2) shape, where n_bonds can vary per frame.

    set read_n_atoms to False when reading old bond reporter outputs which did not contain the n_atoms header.
    """
    data = read_compressed(path)
    if read_n_atoms:
        n_atoms, = struct.unpack("=Q", data[0:8])
        frame_start = 8
    else:
        n_atoms = 0
        frame_start = 0
    frames = []
    while frame_start < len(data):
        n_bonds, = struct.unpack("=Q", data[frame_start:frame_start+8])
        frame_end = frame_start + 8 + 8*n_bonds
        frames.append(
            np.frombuffer(data[frame_start+8:frame_end], dtype=np.uint32)
            .reshape((n_bonds, 2))
        )
        frame_start = frame_end
    n_frames = len(frames)
    return n_frames, n_atoms, frames


def generate_vmd_readable_bonds(
    frames: list[np.ndarray], xtc_path: str, out_path: str
) -> None:
    """
    Given a list of bonds (from a .bonds file, loaded with
    old_formats.read_bonds), and a matching trajectory .xtc file
    (input) write a file in a format that is
    suitable for visualization by daemon.tcl.

    May take a while, reports progress to stdout.
    """

    xtc = molly.XTCReader(xtc_path)
    print("Generating file for vmd bond visualization.")

    comp_obj = zlib.compressobj(6)
    with open(out_path, "wb") as handle:
        for i, frame in enumerate(frames):
            if i % 8 == 0:
                print(f"Processing frame {i} of {len(frames)}")
            xtc_frame = xtc.pop_frame()
            pos = xtc_frame.positions
            box = np.array(
                [
                    xtc_frame.box[0][0],
                    xtc_frame.box[1][1],
                    xtc_frame.box[2][2]
                ],
                dtype=np.float32
            )
            n_atoms = len(pos)
            bonds = [[] for _ in range(n_atoms)]
            crosses = np.any(
                np.abs(
                    pos[frame[:, 0]] - pos[frame[:, 1]]
                ) / box > 0.5,
                axis=-1
            )
            for [i, j], cross in zip(frame, crosses):
                if cross:
                    continue
                bonds[i].append(j)
                bonds[j].append(i)
            handle.write(
                comp_obj.compress(
                    (
                        "{" + "} {".join([
                            " ".join([
                                f"{other_atom}"
                                for other_atom in other_atom_list
                            ]) for other_atom_list in bonds
                        ]) + "} "
                    ).encode("utf-8")
                )
            )
        handle.write(comp_obj.flush())
    print("Generating file finished.")

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
    Parses the old (pre 0.2) .frags format,
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