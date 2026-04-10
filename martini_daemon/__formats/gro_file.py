# Simple gromacs gro file read/write
import numpy as np

from ..__rust import PeriodicBox


def read_gro(path):
    """Reads .gro file at path, returns box, pos, vel.
    box is a PeriodicBox instance.
    pos and vel are float64 numpy arrays.

    If there are no velocities in the gro file, returns None for vel.
    (if even a single velocity is missing, it returns None).
    """
    with open(path) as file:
        file.readline()  # skip title
        try:
            n_atoms = int(file.readline().strip())
        except ValueError:
            raise ValueError(
                "Error parsing gro file: second line should be the number of atoms."
            )
        pos = np.zeros((n_atoms, 3), dtype=np.float64)
        vel = np.zeros((n_atoms, 3), dtype=np.float64)
        for i in range(n_atoms):
            line = file.readline()
            if len(line) < 44:
                raise ValueError(
                    f"Error parsing gro file, line for atom index {i} (zero indexed) is shorter than 44 characters."
                )
            try:
                pos[i][0] = float(line[20:28].strip())
                pos[i][1] = float(line[28:36].strip())
                pos[i][2] = float(line[36:44].strip())
            except ValueError:
                raise ValueError(
                    f"Error parsing gro file, invalid coordinate / floating point number for atom index {i} (zero indexed)."
                )
            if len(line) < 68:
                vel = None
            if vel is None:
                continue
            try:
                vel[i][0] = float(line[44:52].strip())
                vel[i][1] = float(line[52:60].strip())
                vel[i][2] = float(line[60:68].strip())
            except ValueError:
                raise ValueError(
                    f"Error parsing gro file, invalid velocity / floating point number for atom index {i} (zero indexed)."
                )
        last_line = file.readline().strip().split()
        try:
            box_floats = [float(x) for x in last_line]
        except ValueError:
            raise ValueError(
                "Error parsing gro file, the coordinates line contains non numbers. Is the number of atoms correct?"
            )
        assert len(box_floats) >= 3
        box = PeriodicBox.from_gro(*box_floats)
        return box, pos, vel


def write_gro(
    path, title, atom_names, res_names, res_ids, box: PeriodicBox, pos, vel=None
) -> None:
    """Write .gro file at path.

    :param path: path to .gro file
    :param title: title of .gro file
    :param atom_names: list of atom names
    :param res_names: list of residue names
    :param res_ids: list of residue IDs
    :param box: PeriodicBox instance.
    :param pos: numpy array of positions in nm.
    :param vel: numpy array of velocities in nm/picosecond or None.
    """
    title = title.strip()
    assert "\n" not in title, "Title must not contain newlines"
    assert len(atom_names) == len(pos) == len(res_names) == len(res_ids)
    if vel is not None:
        assert len(pos) == len(vel)
    n_atoms = len(atom_names)
    with open(path, "w") as file:
        file.write(f"{title}\n")
        file.write(f"{n_atoms}\n")
        for i in range(n_atoms):
            cpos = pos[i]
            name, resid, res_name = atom_names[i], res_ids[i], res_names[i]
            index = i + 1
            file.write(f"{resid % 100000:5}{res_name:5}{name:>5}")
            file.write(f"{index % 100000:5}{cpos[0]:8.3f}{cpos[1]:8.3f}{cpos[2]:8.3f}")
            if vel is not None:
                cvel = vel[i]
                file.write(f"{cvel[0]:8.4f}{cvel[1]:8.4f}{cvel[2]:8.4f}")
            file.write("\n")
        file.write(box.to_gro())
