import os
import sys
from argparse import ArgumentParser
from random import random
from warnings import warn

import numpy as np

from .__formats import (
    TopTrajReader,
    TopTrajWriter,
    TrajectoryReader,
    TrajectoryWriter,
    read_checkpoint,
)
from .__reporters import FragmentReporter
from .__rust import BondGraph, build_version
from .__simulation import Simulation
from .extra import logo, logo_small


def print_help(topic: None | str = None) -> int:
    match topic:
        case None | "help":
            print("""Martini Daemon CLI
    
Usage: daemon [SUBCOMMAND] [OPTIONS]
Subcommands:
help - print this help message
version - print version information
info - print information about a .toptraj or .chk file
whole - make a trajectory whole using a .toptraj file
select - select atoms and/or frames for .toptraj file

Use `daemon help [SUBCOMMAND]` for more information about a subcommand.""")
        case "version":
            print("""daemon version

Prints version information, as shown in .log files.""")
        case "info":
            print("""daemon info

Print information about a .toptraj or .chk file.
Takes one positional argument -- the path to the file to provide more info about.""")
        case "whole":
            print("""daemon whole

Make a trajectory whole using a .toptraj file.
Takes the following arguments:
-i input trajectory file
-o output trajectory file
-s toptraj file""")
        case "select":
            print("""daemon select

Select atoms and/or frames for a .toptraj file.
Takes the following arguments:
-i input toptraj file
-o output toptraj file
-a/--atoms start:end:step
-f/--frames start:end:step
""")
        case "dist":
            print("""daemon dist

Print distance, angle or dihedrals from a trajectory, based on fragments and graph information, for making distributions.
Takes the following arguments:
-i, --input input trajectory file
-o, --output output file
-f, --frags '.frags' file with detailed fragment listing enabled
-p, --topology topology .top file with graph information #included
-s, --selection Selection pattern, e.g. "monomer:r-s" for a distance between nodes r and s in graph monomer.
""")
        case _:
            print(f"Unknown subcommand: {topic}.")
            return 1
    return 0


def info(args: list[str]) -> int:
    """Get information about daemon files using the `daemon info` CLI."""
    if len(args) < 1:
        print_help("info")
        return 1
    path = args[0]
    if not os.path.isfile(path):
        print(f"File {path} does not exist.")
        return 2
    base, ext = os.path.splitext(path)
    print(f"path:                     {path}")
    match ext:
        case ".chk" | ".chk2" | ".chk3":
            chk = read_checkpoint(path)
            print(f"step:                     {chk.current_step}")
            print(f"time:                     {chk.time_ps:.1f} ps")
            print(f"frame index:              {chk.trajectory_frame}")
            print(f"n_atoms:                  {chk.n_atoms}")
            print(f"n_reactions:              {chk.reactions_so_far}")
            print(f"box vectors (nm):         {chk.box.to_lattice()}")
            print(f"largest abs(pos) (nm):    {np.max(np.abs(chk.pos))}")
            print(f"largest abs(vel) (nm/ps): {np.max(np.abs(chk.vel))}")
        case ".toptraj":
            r = TopTrajReader(path)
            # frame start offsets in file
            tells = [r.tell()]
            while r.skip_frame():
                tells.append(r.tell())
            del tells[-1]
            print(f"toptraj version:          {r.major_version}.{r.minor_version}")
            print(f"sim title:                {r.title}")
            print("initial molecules:")
            for name, count, atoms_per in r.initial_molecules:
                print(f"  {name} {count} ({atoms_per} atoms per molecule)")
            print(f"n_frames:                 {len(tells)}")
            print(f"n_atoms:                  {r.n_atoms}")
            r.seek(tells[0])
            first_frame = r.read_frame()
            r.seek(tells[-1])
            last_frame = r.read_frame()
            first_step = f"{first_frame.sim_step}" if first_frame is not None else "n/a"
            last_step = f"{last_frame.sim_step}" if last_frame is not None else "n/a"
            first_time = (
                f"{first_frame.sim_time:.0f}" if first_frame is not None else "n/a"
            )
            last_time = (
                f"{last_frame.sim_time:.0f}" if last_frame is not None else "n/a"
            )
            print(f"step:                     {first_step}-{last_step}")
            print(f"time:                     {first_time}-{last_time} ps")
        case _:
            print(f"No info for file extension: {ext}.")
            return 1
    return 0


def whole(args: list[str]) -> int:
    """Make a trajectory whole using the `daemon whole` CLI."""
    parser = ArgumentParser(
        prog="daemon whole",
        description="Make a trajectory whole using a .toptraj file",
        add_help=False,
    )
    parser.add_argument("-i", "--input", required=True, help="input trajectory file")
    parser.add_argument("-o", "--output", required=True, help="output trajectory file")
    parser.add_argument(
        "-s", "--topology", required=True, help="topology (toptraj file)"
    )
    parsed_args = parser.parse_args(args)
    inp = parsed_args.input
    oup = parsed_args.output
    toptraj = parsed_args.topology

    if inp == oup:
        print(f"Error: {inp} and {oup} are the same path.")
        return 1
    if not os.path.isfile(inp):
        print(f"{inp} does not exist, or is not a file.")
        return 1
    if not os.path.isfile(toptraj):
        print(f"{toptraj} does not exist, or is not a file.")
        return 1

    r = TopTrajReader(toptraj)
    bonds = []
    i = 0
    while f := r.read_frame():
        i += 1
        bonds.append(f.bonds)
        sys.stdout.write(f"\033[2K\rReading bonds: {i}")
    print()

    n_frames = len(bonds)
    n_atoms = r.n_atoms
    r.close()
    print(f"Read {n_frames} frames from {toptraj}.")

    r = TrajectoryReader(inp)
    w = TrajectoryWriter(oup)

    print()
    for i in range(n_frames):
        sys.stdout.write(f"\033[2K\rProcessing frame: {i + 1}/{n_frames + 1}")
        frame = r.read_frame()
        assert frame is not None
        step, time, pbc, pos, _ = frame
        if len(pos) != n_atoms:
            print(
                f"Number of atoms in {toptraj} ({n_atoms}) and in {inp} ({len(pos)}) does not match."
            )
            return 1
        pos = np.array(pos, dtype=np.float64)

        graph = BondGraph(n_atoms)
        for a, b in bonds[i]:
            graph.add_bond(a, b)
        graph.make_whole(pbc, pos)

        w.write_frame(step, time, pbc, pos)

    r.close()
    w.close()
    print()
    return 0


def dist(args: list[str]) -> int:
    """Analyze distributions."""
    parser = ArgumentParser(
        prog="daemon dist",
        description="Analyze bond, angle or dihedral distributions",
        add_help=False,
    )
    parser.add_argument("-i", "--input", required=True, help="input trajectory file")
    parser.add_argument("-o", "--output", required=True, help="output file")
    parser.add_argument(
        "-f",
        "--frags",
        required=True,
        help=".frags file with detailed fragment listing",
    )
    parser.add_argument(
        "-p",
        "--topology",
        required=True,
        help="topology .top file with graph information #included",
    )
    parser.add_argument(
        "-s",
        "--selection",
        required=True,
        help='Selection pattern, e.g. "monomer:r-s" for a distance between nodes r and s in graph monomer. T',
    )

    parsed_args = parser.parse_args(args)
    inp = TrajectoryReader(parsed_args.input)
    frag_frames = FragmentReporter.read_fragments(parsed_args.frags)
    graphs = Simulation(parsed_args.topology, sim_name=None).top.graphs
    graph_name = parsed_args.selection.split(":")[0].strip()
    node_names = [x.strip() for x in parsed_args.selection.split(":")[1].split("-")]

    if graph_name not in graphs:
        print(f"{graph_name} not defined in {parsed_args.topology}.")
        return 1
    graph = graphs[graph_name]

    for node in node_names:
        if node not in graph.atom_name_to_index:
            print(f"{node} not part of graph {graph_name} in {parsed_args.topology}.")
            return 1
    nodes = [graph.atom_name_to_index[x] for x in node_names]
    if len(nodes) < 2 or len(nodes) > 4:
        print(
            "Specify 2 (distance), 3 (angle) or 4 (dihedral) nodes to write to output file."
        )

    with open(parsed_args.output, "w") as f:
        f.write(f"# Written by daemon dist with arguments {args}\n")
        n_tot = 0
        for frame in frag_frames:
            sys.stdout.write(
                f"\033[2K\rAnalyzing frame: {frame.frame_index} out of {len(frag_frames)}"
            )
            if frame.fragments is None:
                warn(f"Frame {frame.frame_index} has no fragments.")
                continue
            if frame.time_ps is None:
                warn(f"Frame {frame.frame_index} has no simulation time.")
                continue
            traj_frame = inp.read_frame()
            if traj_frame is None:
                print(
                    f"Trajectory {parsed_args.input} does not have enough frames in comparison with {parsed_args.frags}."
                )
                return 2
            step, time, box, pos, vel = traj_frame
            f.write(
                f"# Frame {frame.frame_index}, Step {frame.step}, Box {box.to_lattice()}\n"
            )
            if step != frame.step or not np.isclose(time, frame.time_ps, atol=0.1):
                print(
                    f"Trajectory {parsed_args.input} and .frags file {parsed_args.frags} are not in sync."
                )
                return 2

            for frag in frame.fragments:
                n_tot += 1
                if frag.name != graph_name:
                    continue
                idx = [frag.atoms[x] for x in nodes]
                if any(x == -1 for x in idx):
                    continue
                match idx:
                    case i1, i2:
                        f.write(f"{box.distance(pos[i1], pos[i2])}\n")
                    case i1, i2, i3:
                        f.write(f"{box.angle(pos[i1], pos[i2], pos[i3])}\n")
                    case i1, i2, i3, i4:
                        f.write(f"{box.dihedral(pos[i1], pos[i2], pos[i3], pos[i4])}\n")
        sys.stdout.write(
            f"\033[2K\rAnalyzed {len(frag_frames)} frames, wrote a total of {n_tot} entries.\n"
        )

    return 0


def parse_slice(slice_: str, n: int) -> tuple[int, int, int]:
    """Parse a slice string into start, end and step slice."""
    assert n > 0, "n_frames must be positive"
    if len(slice_) == 0:
        return 0, n, 1
    toks = slice_.split(":")
    if len(toks) == 1:
        i = int(slice_)
        return i, i + 1, 1
    if 2 <= len(toks) <= 3:
        start = int(toks[0]) if toks[0] != "" else 0
        if start < 0:
            start += n
        end = int(toks[1]) if toks[1] != "" else n
        if end < 0:
            end += n
        step = int(toks[2]) if len(toks) == 3 else 1
        return start, end, step
    raise ValueError("Too many ':' in slice, up to 2 expected.")


def select(args: list[str]) -> int:
    """Select frames or atoms using the `daemon select` CLI."""
    parser = ArgumentParser(
        prog="daemon select",
        description="Select atoms and/or frames for a .toptraj or .frags file",
        add_help=False,
    )
    parser.add_argument("-i", "--input", required=True, help="input toptraj or frags file")
    parser.add_argument("-o", "--output", required=True, help="output toptraj or frags file")
    parser.add_argument(
        "-a", "--atoms", required=False, help="atom selection (start:stop:step), only valid for toptraj files"
    )
    parser.add_argument(
        "-f", "--frames", required=False, help="frame selection (start:stop:step)"
    )
    parsed_args = parser.parse_args(args)
    inp = parsed_args.input
    oup = parsed_args.output

    if inp == oup:
        print(f"Error: {inp} and {oup} are the same path.")
        return 1
    if not os.path.isfile(inp):
        print(f"{inp} does not exist, or is not a file.")
        return 1
    _, ext = os.path.splitext(inp)
    _, oup_ext = os.path.splitext(oup)
    if oup_ext != ext:
        print(f"Extension of input/output not the same. Input: {ext}, output: {oup_ext}.")
        return 1
    match ext:
        case ".toptraj":
            r = TopTrajReader(inp)
            # count frames
            tells = [r.tell()]
            while r.skip_frame():
                tells.append(r.tell())
            del tells[-1]
            r.seek(tells[0])
            atom_start, atom_end, atom_step = parse_slice(parsed_args.atoms or "", r.n_atoms)
            new_n_atoms = len(range(atom_start, atom_end, atom_step))
            frame_start, frame_end, frame_step = parse_slice(
                parsed_args.frames or "", len(tells)
            )

            if atom_start > atom_end or atom_step < 0:
                raise ValueError("Atom selection error - must have a positive step.")
            if atom_start > 0 or atom_step != 1:
                # to fix this, a map of old atom -> new atom indices would need to be made, applied to bonds, and made usable
                # in e.g. the vmd plugin
                raise ValueError(
                    "Currently, only atom selections starting at 0 and with step 1 are supported."
                )
            if frame_start > frame_end or frame_step < 0:
                raise ValueError("Frame selection error - must have a positive step.")

            res_names = [r.res_names[j] for j in range(atom_start, atom_end, atom_step)]
            res_ids = [r.res_ids[j] for j in range(atom_start, atom_end, atom_step)]

            if (sel_start := r.title.find("(select")) > -1:
                new_title = r.title[:sel_start] + "(selected multiple times)"
            else:
                new_title = (
                    r.title
                    + f" (select atoms: {atom_start}:{atom_end}:{atom_step}; frames: {frame_start}:{frame_end}:{frame_step})"
                )

            w = TopTrajWriter(oup, new_title, r.initial_molecules, res_names, res_ids)

            new_i = -1
            for new_i, i in enumerate(range(frame_start, frame_end, frame_step)):
                sys.stdout.write(f"\033[2K\rProcessing frame: {i + 1}/{len(tells) + 1}")
                r.seek(tells[i])
                frame = r.read_frame()
                assert frame is not None
                w.new_frame(new_i, frame.sim_step, frame.sim_time, new_n_atoms)

                names = [frame.names[j] for j in range(atom_start, atom_end, atom_step)]
                types = [frame.atom_types[j] for j in range(atom_start, atom_end, atom_step)]
                charges = [frame.charges[j] for j in range(atom_start, atom_end, atom_step)]
                masses = [frame.masses[j] for j in range(atom_start, atom_end, atom_step)]
                w.write_frame_atoms(names, types, charges, masses)

                sel = set(range(atom_start, atom_end, atom_step))
                bonds = [(a, b) for a, b in frame.bonds if a in sel and b in sel]
                w.write_frame_bonds(bonds)
                w.write_frame()
            print(f"\033[2K\rProcessed a total of {new_i + 1} frames.")

        case ".frags":
            if parsed_args.atoms is not None:
                warn("Warning: Argument --atoms is ignored on .frags select.")
            with open(inp) as f_in:
                # sadly we need to know the number of frames in advance
                # luckily due to disk caching this isn't a complete waste
                # to read like this
                n_frames = 0
                while line := f_in.readline():
                    if line == "":
                        break
                    if line.startswith("Step:"):
                        n_frames += 1
                
            with open(oup, "w") as f_out:
                frame_start, frame_end, frame_step = parse_slice(
                    parsed_args.frames or "", n_frames
                )
                last_written = -1
                for frame in FragmentReporter.iter_fragments(inp):
                    if frame.frame_index is None:
                        # detailed mode is off, best guess
                        frame.frame_index = last_written + 1
                    if frame.frame_index < frame_start:
                        continue
                    if frame.frame_index >= frame_end:
                        break
                    last_written += 1
                    if last_written % frame_step != 0:
                        continue
                    f_out.write(
                        frame.serialize()
                    )

        case _:
            print(f"File type {ext} not supported, only .toptraj and .frags is.")
    return 0

def main() -> None:
    """Parse arguments and run the right subcommand of the `daemon` CLI."""
    if len(sys.argv) < 2:
        print_help()
        sys.exit(1)
    subcommand = sys.argv[1]
    args = sys.argv[2:]

    match subcommand:
        case "help" | "-h" | "--help":
            sys.exit(print_help(args[0] if len(args) > 0 else None))
        case "version" | "-v" | "--version":
            print(build_version())
            sys.exit(0)
        case "info":
            sys.exit(info(args))
        case "whole":
            sys.exit(whole(args))
        case "select":
            sys.exit(select(args))
        case "dist":
            sys.exit(dist(args))
        case "logo":
            if random() > 0.5:
                print(logo())
            else:
                print(logo_small())
            sys.exit(0)
        case _:
            print(f"Unknown command: {subcommand}")
            print()
            print_help()
            sys.exit(1)
