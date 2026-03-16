# make trajectories whole, or emit cluster information for mdvwhole

import numpy as np
from multiprocessing import Pool
from tqdm import tqdm
from functools import partial
import molly
from ..__core import make_cluster_frame, make_whole_frame

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
            len(pos),
            pos, box, clusters[frame_i],
        )
        f.positions = pos.flatten()
        writer.write_frame(f)
    writer.close()
