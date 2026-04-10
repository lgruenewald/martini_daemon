from os.path import splitext

from .gro_file import read_gro, write_gro
from .xyz_file import read_xyz, write_xyz


def read_geometry(path):
    """File format implied from extension.

    :param path: path to file
    :return: periodic box (PeriodicBox), positions (numpy float64 array), velocities (numpy float64 array or None).
    """
    prefix, ext = splitext(path)
    match ext:
        case ".gro":
            return read_gro(path)
        case ".xyz":
            return read_xyz(path)
        case _:
            raise ValueError(f"Unknown geometry file format {ext} for {path}.")


def write_geometry(path, title, atom_names, res_names, res_ids, box, pos, vel=None):
    """File format implied from extension.

    :param path: path to file
    :param title: title of geometry
    :param atom_names: list of atoms names
    :param res_ids: list of residues IDs
    :param res_names: list of residues names
    :param box: periodic box
    :param pos: positions
    :param vel: velocities or None
    """
    prefix, ext = splitext(path)
    match ext:
        case ".gro":
            return write_gro(path, title, atom_names, res_names, res_ids, box, pos, vel)
        case ".xyz":
            return write_xyz(path, atom_names, box, pos, vel)
        case _:
            raise ValueError(f"Unknown geometry file format {ext} for {path}.")
