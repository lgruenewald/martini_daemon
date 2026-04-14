# Simple xyz file read/write
import re

import numpy as np
import numpy.typing as npt

from ..__rust import PeriodicBox


def read_xyz(
    path: str,
) -> tuple[PeriodicBox, npt.NDArray[np.float64], npt.NDArray[np.float64] | None]:
    """Reads an extended .xyz file at path, returns box, pos, vel.

    Coordinates in the .xyz file are assumed to be in Angstroms.
    Velocities in the .xyz file are assumed to be in Angstroms/picosecond.

    :param path: path to .xyz file to read.
    :return: The periodic box, positions in nanometers, and velocities in nanometers/picosecond if present.
    """
    # unit conversion constants
    pos_conversion = 0.1
    vel_conversion = 0.1

    with open(path) as file:
        try:
            n_atoms = int(file.readline().strip())
        except ValueError:
            raise ValueError(
                "Error parsing xyz file: first line should be the number of atoms."
            )

        title_line = file.readline().strip()
        match = re.match('Lattice="[^"]*"', title_line)
        if match is None:
            raise ValueError(
                "Error parsing xyz file, the title line contains no Lattice description."
            )
        lattice = title_line[match.start() + 9 : match.end() - 1]
        try:
            box_floats = [float(x) for x in lattice.split()]
        except ValueError:
            raise ValueError(
                "Error parsing xyz file, the Lattice description must be only numbers."
            )
        if len(box_floats) != 9:
            raise ValueError(
                "Error parsing xyz file, the Lattice description must be 9 numbers."
            )

        # 0: x1 1: y1 2: z1 3: x2 4: y2 5: z2 6: x3 7: y3 8: z3
        # we need to remove the useless ones - z2, z1, y1
        del box_floats[5]
        del box_floats[2]
        del box_floats[1]
        box = PeriodicBox.triclinic(*box_floats)

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
                raise ValueError(
                    f"Error parsing xyz file, invalid coordinate / floating point number for atom index {i} (zero indexed)."
                )

            if len(tokens) < 7:
                vel = None
            if vel is not None:
                try:
                    vel[i, 0] = float(tokens[4]) * vel_conversion
                    vel[i, 1] = float(tokens[5]) * vel_conversion
                    vel[i, 2] = float(tokens[6]) * vel_conversion
                except ValueError:
                    raise ValueError(
                        f"Error parsing xyz file, invalid velocity / floating point number for atom index {i} (zero indexed)."
                    )
        return box, pos, vel


def write_xyz(
    path: str,
    atoms: list[str],
    box: PeriodicBox,
    pos: npt.NDArray[np.float64 | np.float32],
    vel: npt.NDArray[np.float64 | np.float32] | None = None,
) -> None:
    """Write an .xyz file.

    Does not round / writes all digits, even ones that are not relevant.

    Writes in Angstrom and Angstroms/picosecond.

    :param path: path to .xyz file to write.
    :param atoms: list of atom names.
    :param box: PeriodicBox, units of nm.
    :param pos: numpy array of positions in nm.
    :param vel: numpy array of velocities in nm/picosecond, or None.
    """
    n_atoms = len(atoms)
    assert len(atoms) == len(pos)
    if vel is not None:
        assert len(pos) == len(vel)
    with open(path, "w") as file:
        file.write(f"{len(atoms)}\n")
        file.write(f'Lattice="{box.to_lattice()}"\n')
        for i in range(n_atoms):
            cpos = pos[i] * 10.0
            name = atoms[i]
            file.write(f"{name} {cpos[0]} {cpos[1]} {cpos[2]}")
            if vel is not None:
                cvel = vel[i] * 10.0
                file.write(f" {cvel[0]} {cvel[1]} {cvel[2]}")
            file.write("\n")
