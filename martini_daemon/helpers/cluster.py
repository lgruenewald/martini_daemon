# make trajectories whole, or emit cluster information for mdvwhole

import numpy as np

def make_clusters(n_atoms: int, bonds_data: list[np.array]):
    """
    Given per-frame bonds data read by bond_reporter.load_bonds, it returns
    per-frame cluster data, in the form of a numpy array of dimension
    (n_frames, n_atoms), type uint32, where every atom is given a numeric
    cluster index. Note: The numbering of clusters across different
    frames is not related.
    """
    n_frames = len(bonds_data)
    res = np.zeros((n_frames, n_atoms), dtype=np.uint32) - 1
    unset = 0xffffffff

    def set_cluster(frame, i, cluster):
        # TODO go through the bond graph and add everything to this cluster
        pass

    for frame in range(n_frames):
        cluster = 0
        for i in range(n_atoms):
            if res[frame][i] != unset:
                # already in a cluster
                continue
            else:
                # no cluster yet
                set_cluster(frame, i, cluster)
                cluster += 1
