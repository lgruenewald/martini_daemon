# Simple gromacs gro file read/write
from .utils import backup_try
import numpy as np


def read_gro(path):
    """
        Reads .gro file at path, returns box, pos, vel.
        box is a python tuple of 3 floating point numbers.
        pos and vel are float64 numpy arrays.

        If there are no velocities in the gro file, returns None
        (if even a single velocity is missing, it returns None).
    """

    with open(path, "r") as file:
        file.readline()  # skip title
        try:
            n_atoms = int(file.readline().strip())
        except ValueError:
            raise ValueError(f"Error parsing gro file: second line should be the number of atoms.")
        pos = np.zeros((n_atoms, 3), dtype=np.float64)
        vel = np.zeros((n_atoms, 3), dtype=np.float64)
        for i in range(n_atoms):
            line = file.readline()
            if len(line) < 44:
                raise ValueError(f"Error parsing gro file, line for atom index {i} (zero indexed) is shorter than 44 characters.")
            try:
                pos[i][0] = float(line[21:28].strip())
                pos[i][1] = float(line[29:36].strip())
                pos[i][2] = float(line[37:44].strip())
            except ValueError:
                raise ValueError(f"Error parsing gro file, invalid coordinate / floating point number for atom index {i} (zero indexed).")
            if len(line) < 68:
                vel = None
            if vel is None:
                continue
            try:
                vel[i][0] = float(line[45:52].strip())
                vel[i][1] = float(line[53:60].strip())
                vel[i][2] = float(line[61:68].strip())
            except ValueError:
                raise ValueError(f"Error parsing gro file, invalid velocity / floating point number for atom index {i} (zero indexed).")
        last_line = file.readline().strip().split()
        try:
            box = list(map(lambda x: float(x), last_line))
        except ValueError:
            raise ValueError(f"Error parsing gro file, the coordinates line contains non numbers. Is the number of atoms correct?")
        assert len(box) >= 3
        for val in box[3:]:
            assert val == 0., "Error parsing gro file, only 90 degree angle pbc are supported."
        return tuple(box[0:3]), pos, vel


def write_gro(path, title, parts, box, pos, vel=None):
    """
        Write .gro file at path, with title, parts (list of tuples
        containing name, resid, resname in this order, everything else
        in the tuple is discarded), box, pos and vel.
    """
    backup_try(path)
    n_atoms = len(parts)
    assert len(parts) == len(pos)
    if vel is not None:
        assert len(pos) == len(vel)
    with open(path, "w") as file:
        file.write(f"{title}\n")
        file.write(f"{len(parts)}\n")
        for i in range(n_atoms):
            cpos = pos[i]
            name, resid, resname, *_ = parts[i]
            index = i + 1
            file.write(f"{resid:5}{resname:5}{name:>5}")
            file.write(f"{index:5}{cpos[0]:8.3f}{cpos[1]:8.3f}{cpos[2]:8.3f}")
            if vel is not None:
                cvel = vel[i]
                file.write(f"{cvel[0]:8.4f}{cvel[1]:8.4f}{cvel[2]:8.4f}")
            file.write("\n")
        file.write(f"{box[0]:.4f} {box[1]:.4f} {box[2]:.4f}\n")
