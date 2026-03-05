# Simple xyz file read/write
from .utils import backup_try
import numpy as np
import re

def read_xyz(path, pos_conversion=0.1, vel_conversion=0.1):
    """
        Reads an extended .xyz file at path, returns box, pos, vel.
        box is a python tuple of 3 floating point numbers.
        pos and vel are float64 numpy arrays.

        If there are no velocities in the xyz file, returns None
        (if even a single velocity is missing, it returns None).

        Coordinates are assumed to be in Angstroms by default.
        Velocities are assumed to be in Angstroms/picosecond by default.

        Custom units can be used by passing pos_conversion and vel_conversion.
        The values in the files are multiplied by these factors.
    """

    with open(path, "r") as file:
        try:
            n_atoms = int(file.readline().strip())
        except ValueError:
            raise ValueError(f"Error parsing xyz file: first line should be the number of atoms.")

        try:
            title_line = file.readline().strip()
            match = re.match('Lattice="[^"]*"', title_line)
            lattice = title_line[match.start():match.end()]
            box = list(map(lambda x: float(x), lattice.split()))
            # 0: x1 1: y1 2: z1 3: x2 4: y2 5: z2 6: x3 7: y3 8: z3
            box = [box[0], box[4], box[8]]
        except:
            raise ValueError(f"Error parsing xyz file, the title line contains no or invalid Lattice description?")
        assert len(box) == 3

        pos = np.zeros((n_atoms, 3), dtype=np.float64)
        vel = np.zeros((n_atoms, 3), dtype=np.float64)
        for i in range(n_atoms):
            line = file.readline()
            tokens = line.split()
            # atomname x y z vx vy vz
            try:
                pos[i, 0] = float(tokens[1]) * pos_conversion
                pos[i, 1] = float(tokens[2]) * pos_conversion
                pos[i, 2] = float(tokens[3]) * pos_conversion
            except ValueError:
                raise ValueError(f"Error parsing xyz file, invalid coordinate / floating point number for atom index {i} (zero indexed).")

            if len(tokens) < 7:
                vel = None
            if vel is not None:
                try:
                    vel[i, 0] = float(tokens[4]) * vel_conversion
                    vel[i, 1] = float(tokens[5]) * vel_conversion
                    vel[i, 2] = float(tokens[6]) * vel_conversion
                except ValueError:
                    raise ValueError(f"Error parsing xyz file, invalid velocity / floating point number for atom index {i} (zero indexed).")
        return tuple(box[0:3]), pos, vel


def write_xyz(path, atoms, box, pos, vel=None):
    """
        Write .xyz file at path, atoms (list of tuples
        containing name as first elem, everything else
        in the tuple is discarded), box, pos and vel.

        Does not round / writes all digits, even ones that
        are not relevant.

        Assumes input box, pos and vel are in nm and nm/picosecond.
        Writes in Angstrom and Angstroms/picosecond.
    """
    backup_try(path)
    n_atoms = len(atoms)
    assert len(atoms) == len(pos)
    if vel is not None:
        assert len(pos) == len(vel)
    with open(path, "w") as file:
        file.write(f"{len(atoms)}\n")
        file.write(f'Lattice="{box[0]*10.} 0.0 0.0 0.0 {box[1]*10.} 0.0 0.0 0.0 {box[2]*10.}\n')
        for i in range(n_atoms):
            cpos = pos[i] * 10.
            name, *_ = atoms[i]
            file.write(f"{name} {cpos[0]} {cpos[1]} {cpos[2]}")
            if vel is not None:
                cvel = vel[i] * 10.
                file.write(f" {cvel[0]} {cvel[1]} {cvel[2]}")
            file.write("\n")
