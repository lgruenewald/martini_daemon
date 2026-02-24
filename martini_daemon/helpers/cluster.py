# make trajectories whole, or emit cluster information for mdvwhole

import numpy as np
from multiprocessing import Pool
from tqdm import tqdm
from functools import partial
import molly


def make_cluster_frame(
    frame, n_atoms
):
    """
    Given a list of bonds (frame), it creates an (n_atoms,) shape uint32
    numpy array of clusters. Cluster 0 = unbonded. All other numbers =
    contiguous bonded graph within that numbered cluster. Cluster numbers
    are arbitrary and non-contiguous.

    See make_clusters for doing it for the whole trajectory, in parallel
    and with tqdm.
    """
    res = np.zeros((n_atoms,), dtype=np.uint32)
    last_cluster = 0
    # iteration 1 - over all bonds
    for bond in frame:
        atom_i = bond[0]
        atom_j = bond[1]
        match (res[atom_i], res[atom_j]):
            case (0, 0):
                last_cluster += 1
                res[atom_i] = last_cluster
                res[atom_j] = last_cluster
            case (ci, 0):
                res[atom_j] = ci
            case (0, cj):
                res[atom_i] = cj
            case (ci, cj):
                # both already have clusters assigned
                # different scenarios possible
                # 1. same cluster => no problem
                if ci == cj:
                    continue
                # 2. diff cluster => let's merge them
                larger = max(ci, cj)
                smaller = min(ci, cj)
                res[res == larger] = smaller
    return res


def make_clusters(
    n_frames, n_atoms, bonds
):
    """
    Given per-frame bonds data read by bond_reporter.load_bonds, it returns
    per-frame cluster data, in the form of a list of numpy arrays of dimension
    (n_atoms,), type uint32, where every atom is given a numeric
    cluster index.

    Arguments:
    n_frames, n_bonds, bonds

    Note: The numbering of clusters across different frames is not related.
    Cluster numbers are not contiguous either.

    Note2: Atoms with no bonds are all in cluster 0, actual clusters start
    at 1.

    Note3: Runs on all cores, uses tqdm.
    """

    with Pool() as p:
        res = [x for x in tqdm(
            p.imap(
                partial(make_cluster_frame, n_atoms=n_atoms),
                bonds
            ),
            desc="clustering",
            total=n_frames
        )]

    return np.vstack(res)


def make_whole_frame(
    pos, box, clus, n_atoms
):
    """
    Makes a single frame whole, based on cluster info. Cluster 0 is ignored.
    Modifies pos in place.

    Args:
    pos, box, clus, n_atoms
    """
    root = {}
    for atom_i in range(n_atoms):
        ci = clus[atom_i]
        if ci == 0:
            continue
        # lowest indexed atom will become "root"
        if root.get(ci) is None:
            root[ci] = pos[atom_i]
        # everything else will be moved by the amount of pbc's required
        # to bring it closest to "root"
        else:
            pos[atom_i] += np.floor(
                (root.get(ci) - pos[atom_i]) / box + 0.5
            ) * box


def make_whole(
    clusters,
    from_xtc: str,
    to_xtc: str
):
    """
    Based on clusters generated from make_clusters, it reads from_xtc
    and writes to_xtc, where the clusters have been made whole.
    """
    reader = molly.XTCReader(from_xtc)
    writer = molly.XTCWriter(to_xtc)

    n_frames = len(clusters)

    for frame_i in tqdm(range(n_frames), desc="making whole"):
        f = reader.pop_frame()
        pos = f.positions
        box = np.array([
            f.box[0][0],
            f.box[1][1],
            f.box[2][2]
        ])
        make_whole_frame(
            pos, box, clusters[frame_i],
            len(pos)
        )
        f.positions = pos.flatten()
        writer.write_frame(f)
    writer.close()
