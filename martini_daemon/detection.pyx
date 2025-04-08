# vectorized detection algorithm
from .utils cimport pdist, pdihedral, pcos_angle, i64
from .fragment import Fragment
import numpy as np
from freud.box import Box
from freud.locality import NeighborList, AABBQuery
import random

cdef i64 get_atom(frag, atom_name) except -2:
    if atom_name in frag.particles:
        return frag.particles[atom_name]
    elif atom_name in frag.opt:
        if frag.opt[atom_name] is None:
            return -1
        else:
            return frag.opt[atom_name]
    else:
        raise ValueError(f"Invalid atom name {atom_name} for frag {frag.name}")

cdef i64 index_pair(frags, i64 frag_index, atom_name) except -2:
    return get_atom(frags[frag_index], atom_name)

cdef bint detection_fragments(
    top, 
    reactants: list[Fragment],
    rx,
    double[:, :] pos,
    double[:] box
):
    """Returns True if frag1 and frag2 fulfill constraints specified in rx

    Returns False if they should not react
    """
    # the order of checks should be from fastest to slowest to maximize
    # performance
    # simple rules check
    # fragments that were removed cannot react any more

    # overlapping fragments can never react:
    cdef set[int] pset = set()
    # TODO further cythonize
    for r in reactants:
        for p in r.particles.values():
            if p in pset:
                return False
            pset.add(p)
        for p in filter(lambda x: x is not None, r.opt.values()):
            if p in pset:
                return False
            pset.add(p)

    # limiter checks, first for performance
    if rx.global_limit is not None and rx.global_counter >= rx.global_limit:
        return False
    if random.random() > rx.probability:
        return False

    cdef double dist
    cdef long idi
    cdef long idj
    cdef double rmax
    # position dependent checks
    # writing it in this style to avoid using iterators
    for i in range(len(rx.distance_max)):
        (idi, namei), (idj, namej), rmax = rx.distance_max[i]
        init1 = index_pair(reactants, idi, namei)
        init2 = index_pair(reactants, idj, namej)
        if init1 == -1 or init2 == -1:
            # missing optional atoms
            continue
        dist = pdist(pos[init1], pos[init2], box)
        if dist > rmax:
            return False

    cdef double rmin
    for i in range(len(rx.distance_min)):
        (idi, namei), (idj, namej), rmin = rx.distance_min[i]
        init1 = index_pair(reactants, idi, namei)
        init2 = index_pair(reactants, idj, namej)
        if init1 == -1 or init2 == -1:
            # missing optional atoms
            continue
        dist = pdist(pos[init1], pos[init2], box)
        if dist < rmin:
            return False

    cdef long idk
    cdef double cos_min
    cdef double cos_max
    for i in range(len(rx.angle_limits)):
        (idi, namei), (idj, namej), (idk, namek), cos_min, cos_max = rx.angle_limits[i]
        p1 = index_pair(reactants, idi, namei)
        p2 = index_pair(reactants, idj, namej)
        p3 = index_pair(reactants, idk, namek)
        if p1 == -1 or p2 == -1 or p3 == -1:
            # missing optional atoms
            continue
        # the particle positions of particle i, j, k
        cos = pcos_angle(pos[p1], pos[p2], pos[p3], box)
        if cos <= cos_min and cos >= cos_max:
            # inverted comparison because cosine is a constantly decreasing
            # function, cos_min is the minimum angle => max cosine value
            # cos_max is the maximum angle => min cosine value
            return False

    cdef long idl
    cdef double min
    cdef double max
    for i in range(len(rx.dihedral_limits)):
        (idi, namei), (idj, namej), (idk, namek), (idl, namel), min, max = rx.dihedral_limits[i]
        # particle positions
        p1 = index_pair(reactants, idi, namei)
        p2 = index_pair(reactants, idj, namej)
        p3 = index_pair(reactants, idk, namek)
        p4 = index_pair(reactants, idl, namel)
        if p1 == -1 or p2 == -1 or p3 == -1 or p4 == -1:
            # missing optional atoms
            continue
        theta = pdihedral(pos[p1], pos[p2], pos[p3], pos[p4], box)
        if theta >= min and theta <= max:
            return False

    rx.global_counter += 1
    return True


cdef update_query(top, double[:] box, double[:, :] pos):
    # TODO cythonize
    if len(top.frag_list) == 0:
        return
    pos_filtered = []
    nlist_to_frag = {}
    for id, frag in top.frag_list.items():
        filtered_id = len(pos_filtered)
        nlist_to_frag[filtered_id] = id
        atom_map_entry = top.neighbor_atom_map.get(frag.name)
        if atom_map_entry is not None:
            name, _ = atom_map_entry
            part_id = get_atom(frag, name)
            if part_id == -1:
                raise ValueError(f"Internal error: invalid r_max nlist atom for {frag.name}: {name}.")
            pos_filtered.append(pos[part_id])
        else:
            pos_filtered.append(pos[frag.index_atom(0)])
    pos_filtered = np.array(pos_filtered)
    query_box = Box(box[0], box[1], box[2])
    query = NeighborList()
    query = AABBQuery(query_box, pos_filtered)
    top.query = query
    top.nlist_to_frag = nlist_to_frag

cpdef detection(
top,
int i,
double[:] box,
double[:, :] pos,
):
    top.pre_detection(i)
    update_query(top, box, pos)

    reactions = []
    skip: set[int] = set()

    for i, frag_i in top.frag_list.items():
        if i in skip:
            continue
        uni_rx = top.reactions.get((frag_i.name, None, None, None))
        if uni_rx is not None and len(uni_rx) > 0:
            for rx in uni_rx:
                if top.detection([frag_i], rx, pos, box):
                    skip.add(i)
                    reactions.append(([frag_i], rx))
                    break
            if i in skip:
                continue
        if not top.r_continue.get((frag_i.name, None, None)):
            continue
        name, _ = top.neighbor_atom_map[frag_i.name]
        part_id = get_atom(frag_i, name)
        if part_id == -1:
            raise ValueError("Internal error: not found in detection_all")
        pos_i = np.array([pos[part_id]])
        nlist = top.query.query(
            pos_i, {"r_max": top.nlist_cutoff}
        ).toNeighborList()
        for _, neigh_j in nlist[:]:
            j = top.nlist_to_frag[neigh_j]
            # i can become skipped during a nested call
            if i in skip:
                break
            if j in skip or i == j:
                continue
            frag_j = top.frag_list[j]
            if frag_i.name == frag_j.name and j < i:
                # don't double count same-type reactions
                continue
            bi_rx = top.reactions.get((frag_i.name, frag_j.name, None, None))
            if bi_rx is not None and len(bi_rx) > 0:
                for rx in bi_rx:
                    if top.detection([frag_i, frag_j], rx, pos, box):
                        skip.add(i)
                        skip.add(j)
                        reactions.append(([frag_i, frag_j], rx))
                        break
            if not top.r_continue.get((frag_i.name, frag_j.name, None)):
                continue
            for _, neigh_k in nlist[:]:
                k = top.nlist_to_frag[neigh_k]
                if i in skip or j in skip:
                    break
                if k in skip or i == k or j == k:
                    continue
                frag_k = top.frag_list[k]
                if frag_i.name == frag_k.name and k < i:
                    continue
                if frag_j.name == frag_k.name and k < j:
                    continue
                tri_rx = top.reactions.get((frag_i.name, frag_j.name, frag_k.name, None))
                if tri_rx is not None and len(tri_rx) > 0:
                    for rx in tri_rx:
                        if top.detection([frag_i, frag_j, frag_k], rx, pos, box):
                            skip.add(i)
                            skip.add(j)
                            skip.add(k)
                            reactions.append(([frag_i, frag_j, frag_k], rx))
                            break
                if not top.r_continue.get((frag_i.name, frag_j.name, frag_k.name)):
                    continue
                for _, neigh_l in nlist[:]:
                    l = top.nlist_to_frag[neigh_l]
                    if i in skip or j in skip or k in skip:
                        break
                    if l in skip or i == l or j == l or k == l:
                        continue
                    frag_l = top.frag_list[l]
                    if frag_i.name == frag_l.name and l < i:
                        continue
                    if frag_j.name == frag_l.name and l < j:
                        continue
                    if frag_k.name == frag_l.name and l < k:
                        continue
                    tetra_rx = top.reactions.get((frag_i.name, frag_j.name, frag_k.name, frag_l.name))
                    if tetra_rx is not None and len(tetra_rx) > 0:
                        for rx in tetra_rx:
                            if top.detection([frag_i, frag_j, frag_k, frag_l], rx, pos, box):
                                skip.add(i)
                                skip.add(j)
                                skip.add(k)
                                skip.add(l)
                                reactions.append(([frag_i, frag_j, frag_k, frag_l], rx))
                                break

    return reactions

