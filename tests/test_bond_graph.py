import numpy as np
from martini_daemon import PeriodicBox, BondGraph
import random
import pytest
from time import time

def test_make_whole():
    n_atoms = 50000
    n_mol = 100
    atoms_per_mol = n_atoms // n_mol
    bond_length = 0.5

    box = PeriodicBox.orthogonal(5., 6., 7.)
    box_np = np.array([box.a[0], box.b[1], box.c[2]])
    bg = BondGraph(n_atoms)
    pos = np.empty((n_atoms, 3), dtype=np.float64)

    for mol_idx in range(n_mol):
        for atom_in_mol in range(atoms_per_mol):
            atom_idx = mol_idx * atoms_per_mol + atom_in_mol

            if atom_in_mol == 0:
                pos[atom_idx] = np.random.random(3) * box_np
            else:
                disp = np.random.random(3)
                disp *= bond_length / np.linalg.norm(disp)
                pos[atom_idx] = pos[atom_idx - 1] + disp
                bg.add_bond(atom_idx-1, atom_idx)
                if atom_in_mol >= 2 and np.random.random() > 0.5:
                    bg.add_bond(atom_idx-2, atom_idx)
                if atom_in_mol >= 3 and np.random.random() > 0.5:
                    bg.add_bond(atom_idx-3, atom_idx)
                if atom_in_mol >= 4 and np.random.random() > 0.5:
                    bg.add_bond(atom_idx-4, atom_idx)

    pos_within = pos.copy()
    box.move_all_within(pos_within)
    pos_whole = pos_within.copy()
    bg.make_whole(box, pos_whole)

    for (i, j) in bg.to_list():
        assert not box.crosses_box(pos_whole[i], pos_whole[j]), (
            f"Bond {i} to {j}\n pos {pos[i]}, {pos[j]};\n pos_within {pos_within[i]}, {pos_within[j]};\n pos_whole {pos_whole[i]}, {pos_whole[j]};"
        )

    assert np.allclose(pos_whole, pos)

@pytest.mark.parametrize("bonds_per_mol", [250, 500, 1500])
def test_reachable_from(bonds_per_mol):
    n_atoms = 50000
    min_n_mol = 100
    atoms_per_mol = n_atoms // min_n_mol
    not_taken = list(range(n_atoms))

    # some of these molecules may end up ultimately unconnected
    molecules = []

    for i in range(min_n_mol):
        molecules.append([])
        for j in range(atoms_per_mol):
            c = random.choice(not_taken)
            molecules[i].append(c)
            not_taken.remove(c)

    start = time()
    bg = BondGraph(n_atoms)
    for mol_idx in range(min_n_mol):
        for _ in range(bonds_per_mol):
            i = random.choice(molecules[mol_idx])
            while (j := random.choice(molecules[mol_idx])) == i:
                pass
            bg.add_bond(
                i, j
            )
    print(f"bg construction took {time() - start} s.")

    # old cluster code
    start = time()
    clus = np.zeros((n_atoms,), dtype=np.uint32)
    last_cluster = 0
    # iteration 1 - over all bonds
    for bond in bg.to_list():
        atom_i = bond[0]
        atom_j = bond[1]
        match (clus[atom_i], clus[atom_j]):
            case (0, 0):
                last_cluster += 1
                clus[atom_i] = last_cluster
                clus[atom_j] = last_cluster
            case (ci, 0):
                clus[atom_j] = ci
            case (0, cj):
                clus[atom_i] = cj
            case (ci, cj):
                # both already have clusters assigned
                # different scenarios possible
                # 1. same cluster => no problem
                if ci == cj:
                    continue
                # 2. diff cluster => let's merge them
                larger = max(ci, cj)
                smaller = min(ci, cj)
                clus[clus == larger] = smaller
    print(f"clustering took {time() - start} s.")

    # for each atom, get the whole molecule it belongs to and compare with clustering
    start = time()
    visited = set()
    for i in range(n_atoms):
        if i in visited:
            continue
        mol = bg.reachable_from({i})
        visited |= mol
        mol_clus = np.nonzero(clus == clus[i])[0]
        if clus[i] == 0:
            assert len(mol) == 1
        else:
            assert len(mol) == len(mol_clus), f"mol {mol} clus {mol_clus}"
            assert set(mol) == set(mol_clus)
    print(f"comparison took {time() - start} s.")
