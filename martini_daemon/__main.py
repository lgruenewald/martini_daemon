from argparse import ArgumentParser
import sys
import os
import numpy as np

from .__rust import build_version, BondGraph
from .__formats import read_checkpoint, TopTrajReader, TrajectoryReader, TrajectoryWriter, TopTrajWriter

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
        case _:
            print(f"Unknown subcommand: {topic}.")
            return 1
    return 0

def info(args: list[str]) -> int:
    """The `daemon info` CLI."""
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
        case ".chk":
            chk = read_checkpoint(path)
            print(f"step:                     {chk.current_step}")
            print(f"time:                     {chk.time_ps} ps")
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
            print(f"step:                     {first_frame.sim_step}-{last_frame.sim_step}")
            print(f"time:                     {first_frame.sim_time:.0f}-{last_frame.sim_time:.0f} ps")
        case _:
            print(f"No info for file extension: {ext}.")
            return 1
    return 0

def whole(args: list[str]) -> int:
    """The `daemon whole` CLI."""
    parser = ArgumentParser(
        prog="daemon whole",
        description="Make a trajectory whole using a .toptraj file",
        add_help=False
    )
    parser.add_argument("-i", "--input", required=True, help="input trajectory file")
    parser.add_argument("-o", "--output", required=True, help="output trajectory file")
    parser.add_argument("-s", "--topology", required=True, help="topology (toptraj file)")
    parsed_args = parser.parse_args(args)
    inp = parsed_args.i
    oup = parsed_args.o
    toptraj = parsed_args.s

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
            print(f"Number of atoms in {toptraj} ({n_atoms}) and in {inp} ({len(pos)}) does not match.")
            return 1
        pos = np.array(pos, dtype=np.float64)

        graph = BondGraph(n_atoms)
        for a, b in bonds[i]:
            graph.add_bond(a, b)
        graph.make_whole(pbc, pos)

        w.write_frame(
            step, time, pbc, pos
        )

    r.close()
    w.close()
    print()
    return 0

def parse_slice(slice_: str, n: int) -> tuple[int, int, int]:
    """Parse a slice string into start, end and step slice."""

    if len(slice_) == 0:
        return 0, n, 1
    else:
        toks = slice_.split(":")
        if len(toks) == 1:
            i = int(slice_)
            return i, i+1, 1
        elif 2 <= len(toks) <= 3:
            start = int(toks[0]) if toks[0] != "" else 0
            if start < 0:
                start += n
            end = int(toks[1]) if toks[1] != "" else n
            if end < 0:
                end += n
            step = int(toks[2]) if len(toks) == 3 else 1
            return start, end, step
        else:
            raise ValueError("Too many ':' in slice, up to 2 expected.")

def select(args: list[str]) -> int:
    """The `daemon select` CLI."""
    parser = ArgumentParser(
        prog="daemon select",
        description="Select atoms and/or frames for a .toptraj file",
        add_help=False
    )
    parser.add_argument("-i", "--input", required=True, help="input toptraj file")
    parser.add_argument("-o", "--output", required=True, help="output toptraj file")
    parser.add_argument("-a", "--atoms", required=False, help="atom selection (start:stop:step)")
    parser.add_argument("-f", "--frames", required=False, help="frame selection (start:stop:step)")
    parsed_args = parser.parse_args(args)
    inp = parsed_args.i
    oup = parsed_args.o

    if inp == oup:
        print(f"Error: {inp} and {oup} are the same path.")
        return 1
    if not os.path.isfile(inp):
        print(f"{inp} does not exist, or is not a file.")
        return 1

    r = TopTrajReader(inp)
    # count frames
    tells = [r.tell()]
    while r.skip_frame():
        tells.append(r.tell())
    del tells[-1]
    r.seek(tells[0])

    atom_start, atom_end, atom_step = parse_slice(parsed_args.atoms or "", r.n_atoms)
    new_n_atoms = len(range(atom_start, atom_end, atom_step))
    frame_start, frame_end, frame_step = parse_slice(parsed_args.frames or "", len(tells))

    if atom_start > atom_end or atom_step < 0:
        print("Atom selection error - must have a positive step.")
    if frame_start > frame_end or frame_step < 0:
        print("Frame selection error - must have a positive step.")

    res_names = [
        r.res_names[j]
        for j in range(atom_start, atom_end, atom_step)
    ]
    res_ids = [
        r.res_ids[j]
        for j in range(atom_start, atom_end, atom_step)
    ]


    if (sel_start := r.title.find("(select")) > -1:
        new_title = r.title[:sel_start] + "(selected multiple times)"
    else:
        new_title = r.title + f" (select atoms: {atom_start}:{atom_end}:{atom_step}; frames: {frame_start}:{frame_end}:{frame_step})"

    w = TopTrajWriter(
        oup,
        new_title,
        r.initial_molecules,
        res_names,
        res_ids
    )

    new_i = -1
    for new_i, i in enumerate(range(frame_start, frame_end, frame_step)):
        sys.stdout.write(f"\033[2K\rProcessing frame: {i + 1}/{len(tells) + 1}")
        r.seek(tells[i])
        frame = r.read_frame()
        assert frame is not None
        w.new_frame(
            new_i,
            frame.sim_step,
            frame.sim_time,
            new_n_atoms
        )

        names = [
            frame.names[j]
            for j in range(atom_start, atom_end, atom_step)
        ]
        types = [
            frame.atom_types[j]
            for j in range(atom_start, atom_end, atom_step)
        ]
        charges = [
            frame.charges[j]
            for j in range(atom_start, atom_end, atom_step)
        ]
        masses = [
            frame.masses[j]
            for j in range(atom_start, atom_end, atom_step)
        ]
        w.write_frame_atoms(
            names, types, charges, masses
        )

        sel = set(range(atom_start, atom_end, atom_step))
        bonds = [
            (a, b)
            for a, b in frame.bonds
            if a in sel and b in sel
        ]
        w.write_frame_bonds(
            bonds
        )
        w.write_frame()
    print(f"\033[2K\rProcessed a total of {new_i+1} frames.")

    return 0

def logo():
    print(r"""@                                                         %                     
@@                  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ %                      
@ @        @@@@@@                                       %   @@@@               @
   @       @@@@                                        %   @@@@@             @@@
 @   @      @@          @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ @@  @            @@ @@
 @     @@@@@@@@   ......... ...................}......}     @@          @@@   @@
  @       @@@@ @ ..... .......................]........}...@@@@@@@@@@@@       @ 
   @      @ @@  @@  . .........................}.......]..@ @@               @@ 
    @     @  @@   @ ........................]<<.:<___<].)@ @@@              @@  
     @     @  @@   @ . ....................)......]....@   @@             @@    
      @    @  @@    @  ....................].......[.=@    @@            @@     
        @@  @@@@     @@.....................}......^@@     @@          @@       
           @@@@        @....................<{)_{~@@        @@@@@@@@@@@         
                         @.................%~... @                              
                          @...............%.. ..@                               
                           @.............%....@@                                 
                             @..........%..@@@                                  
                               @@@@@@@@@@@ @@                                   
                                @         @                                     
                                 @@     @@                                      
                                   @@@@@                                        
                                   @  @@                                        
                                   @  @@                                        
                                   @  @@                                        
                                   @  @@                                        
                                   @  @@                                        
                                   @  @@                                        
                                   @  @@                                        
                                   @  @@                                @@@     
                                   @  @@                              @@@@      
                                   @  @@              @@          @@@@@@@       
                                   @  @@         @@@@@@@@@@  @@@@@@@@@@@        
                                   @  @@     @@@@@@      @@      @@@@@@         
                                   @  @@  @@@@@@@         @@   @@@@@@@@         
                                   @  @@@@@@@@             @@@@@   @@@          
                               @@@@@  @@@@@@@                       @           
                         @@@@     @     @@     @@@                              
                      @@        @        @@@       @                            
                     @@        @@@@@@@@@@@@         @                           
                      @@@                           @                           
                        @@@@@@@                @@@@@                            
                              @@@@@@@@@@@@@@@@@@                                
                                                                                
    __  __            _   _       _    ____                                     
   |  \/  | __ _ _ __| |_(_)_ __ (_)  |  _ \  __ _  ___ _ __ ___   ___  _ __    
   | |\/| |/ _` | '__| __| | '_ \| |  | | | |/ _` |/ _ \ '_ ` _ \ / _ \| '_ \   
   | |  | | (_| | |  | |_| | | | | |  | |_| | (_| |  __/ | | | | | (_) | | | |  
   |_|  |_|\__,_|_|   \__|_|_| |_|_|  |____/ \__,_|\___|_| |_| |_|\___/|_| |_|  
""")

def main():
    """The `daemon` CLI main function."""

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
        case "logo":
            logo()
            sys.exit(0)
        case _:
            print(f"Unknown command: {subcommand}")
            print()
            print_help()
            sys.exit(1)
