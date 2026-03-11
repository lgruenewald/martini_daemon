from .system import System
from .bonded_force import BondedForce
import numpy as np
from ..__rust import PeriodicBox

def collect_bonds(system: System, filters: list[str]):
    n = system.atom_count()
    bonds = []
    for force in system.get_forces():
        if not issubclass(type(force), BondedForce):
            continue
        # do we include this force
        do_force = False
        if filters is None:
            do_force = True
        else:
            for filt in filters:
                if force.passes_filter(filt):
                    do_force = True
                    break
        for members, _ in force.iterate_bonds():
            if force.passes_filter("vsite"):
                # hardcoded special case, modeled as vsite bonded to all constructing particles
                i = members[0]
                for j in members[1:]:
                    if i == j:
                        continue
                    bonds.append((i, j))
            else:
                # modeled as each particle bonded to the next one
                for i, j in zip(members[:-1], members[1:]):
                    if i == j:
                        continue
                    bonds.append((i, j))
    return bonds

def collect_bonds_for_whole(system: System):
    n = system.atom_count()
    bonds = []
    for force in system.get_forces():
        # do we include this force
        if not issubclass(type(force), BondedForce) or force.uses_pbc():
            continue
        for members, _ in force.iterate_bonds():
            i = members[0]
            for j in members[1:]:
                if i == j:
                    continue
                bonds.append((i, j))
    return bonds

def make_cluster_frame(n_atoms, frame):
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


def make_whole_frame(n_atoms, pos, box: PeriodicBox, clus):
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
            pos[atom_i] = box.move_to(root[ci], pos[atom_i])
