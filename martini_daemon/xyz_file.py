# Simple xyz file read/write
from .utils import backup_try
import numpy as np


def read_xyz(path):
    """
        Reads .xyz file at path, returns box, pos, vel.
        box is a python tuple of 3 floating point numbers.
        pos and vel are float64 numpy arrays.

        If there are no velocities in the xyz file, returns None
        (if even a single velocity is missing, it returns None).
    """

    # TODO fix units, it should be angstrom<=>nm converting
    with open(path, "r") as file:
        try:
            n_atoms = int(file.readline().strip())
        except ValueError:
            raise ValueError(f"Error parsing xyz file: first line should be the number of atoms.")

        title_line = file.readline().strip().split()
        try:
            box = list(map(lambda x: float(x), title_line))
        except ValueError:
            raise ValueError(f"Error parsing xyz file, the title line contains no pbc description?")
        assert len(box) >= 3
        for val in box[3:]:
            assert val == 0., "Error parsing xyz file, only 90 degree angle pbc are supported."

        pos = np.zeros((n_atoms, 3), dtype=np.float64)
        vel = np.zeros((n_atoms, 3), dtype=np.float64)
        for i in range(n_atoms):
            line = file.readline()
            tokens = line.split()
            # atomname x y z vx vy vz
            try:
                pos[i, 0] = float(tokens[1])
                pos[i, 1] = float(tokens[2])
                pos[i, 2] = float(tokens[3])
            except ValueError:
                raise ValueError(f"Error parsing xyz file, invalid coordinate / floating point number for atom index {i} (zero indexed).")

            if len(tokens) < 7:
                vel = None
            if vel is not None:
                try:
                    vel[i, 0] = float(tokens[4])
                    vel[i, 1] = float(tokens[5])
                    vel[i, 2] = float(tokens[6])
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
    """
    backup_try(path)
    n_atoms = len(atoms)
    assert len(atoms) == len(pos)
    if vel is not None:
        assert len(pos) == len(vel)
    with open(path, "w") as file:
        file.write(f"{len(atoms)}\n")
        file.write(f"{box[0]} {box[1]} {box[2]}\n")
        for i in range(n_atoms):
            cpos = pos[i]
            name, *_ = atoms[i]
            file.write(f"{name} {cpos[0]} {cpos[1]} {cpos[2]}")
            if vel is not None:
                cvel = vel[i]
                file.write(f" {cvel[0]} {cvel[1]} {cvel[2]}")
            file.write("\n")
