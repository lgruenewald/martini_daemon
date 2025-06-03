# vectorized detection algorithm
from .utils cimport pdist, pdihedral, pcos_angle, i64
from .fragment import Fragment
import numpy as np
from freud.box import Box
from freud.locality import NeighborList, AABBQuery
import random
import cython

@cython.boundscheck(False)
@cython.wraparound(False)
cdef bint detection_one(
    reactants: list[Fragment],
    rx,
    double[:, :] pos,
    double[:] box,
    double absolute_rate
):
    """Returns True if frag1 and frag2 fulfill constraints specified in rx

    Returns False if they should not react
    """
    # the order of checks should be from fastest to slowest to maximize
    # performance
    # simple rules check
    # fragments that were removed cannot react any more

    # overlapping fragments can never react:
    # TODO skiplists?
    cdef set[i64] pset = set()
    # TODO further cythonize
    for r in reactants:
        for p in r.atoms:
            if p == -1:
                continue
            if p in pset:
                return False
            pset.add(p)

    # all the variables we actually want typed
    cdef double dist
    cdef double min
    cdef double max
    cdef double theta
    cdef double cos_theta
    cdef double prob
    cdef i64 p1
    cdef i64 p2
    cdef i64 p3
    cdef i64 p4

    # position dependent checks
    # writing it in this style to avoid using iterators
    for i in range(len(rx.distance_max)):
        idi, atomi, idj, atomj, max = rx.distance_max[i]
        p1 = reactants[idi].atoms[atomi]
        p2 = reactants[idj].atoms[atomj]
        if p1 == -1 or p2 == -1:
            # missing optional atoms
            continue
        dist = pdist(pos[p1], pos[p2], box)
        if dist > max:
            return False

    for i in range(len(rx.distance_min)):
        idi, atomi, idj, atomj, min = rx.distance_min[i]
        p1 = reactants[idi].atoms[atomi]
        p2 = reactants[idj].atoms[atomj]
        if p1 == -1 or p2 == -1:
            # missing optional atoms
            continue
        dist = pdist(pos[p1], pos[p2], box)
        if dist < min:
            return False

    for i in range(len(rx.angle_limits)):
        idi, atomi, idj, atomj, idk, atomk, min, max = rx.angle_limits[i]
        p1 = reactants[idi].atoms[atomi]
        p2 = reactants[idj].atoms[atomj]
        p3 = reactants[idk].atoms[atomk]
        if p1 == -1 or p2 == -1 or p3 == -1:
            # missing optional atoms
            continue
        # the particle positions of particle i, j, k
        cos_theta = pcos_angle(pos[p1], pos[p2], pos[p3], box)
        if cos_theta <= min and cos_theta >= max:
            # inverted comparison because cosine is a constantly decreasing
            # function, cos_min is the minimum angle => max cosine value
            # cos_max is the maximum angle => min cosine value
            return False

    for i in range(len(rx.dihedral_limits)):
        idi, atomi, idj, atomj, idk, atomk, idl, atoml, min, max = rx.dihedral_limits[i]
        # particle positions
        p1 = reactants[idi].atoms[atomi]
        p2 = reactants[idj].atoms[atomj]
        p3 = reactants[idk].atoms[atomk]
        p4 = reactants[idl].atoms[atoml]
        if p1 == -1 or p2 == -1 or p3 == -1 or p4 == -1:
            # missing optional atoms
            continue
        theta = pdihedral(pos[p1], pos[p2], pos[p3], pos[p4], box)
        if theta >= min and theta <= max:
            return False

    # rate limiting, as a last condition
    if rx.relative_rate is not None:
        rx.reaction_counter += 1
        if absolute_rate <= 0.:
            # negative stands for uninitialized / warmup phase globally
            return False
        if rx.observed_rate is None:
            # warmup phase for this reaction
            return False
        prob = absolute_rate / rx.observed_rate
        assert prob >= 0., "internal error, probability below 0"
        assert prob <= 1., "internal error, probability above 1"
        if prob < random.random():
            # prob=0.0 - always false
            # prob=1.0 - always true
            return False
    
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
            atomi, _ = atom_map_entry
            part_id = frag.atoms[atomi]
            if part_id == -1:
                raise ValueError(f"Internal error: invalid r_max nlist atom for {frag.name}: {atomi}.")
            pos_filtered.append(pos[part_id])
        else:
            pos_filtered.append(pos[frag.index_atom(0)])
    pos_filtered = np.array(pos_filtered)
    query_box = Box(box[0], box[1], box[2])
    query = NeighborList()
    query = AABBQuery(query_box, pos_filtered)
    top.query = query
    top.nlist_to_frag = nlist_to_frag

cpdef detection_all(
top,
double[:] box,
double[:, :] pos,
):
    update_query(top, box, pos)

    reactions = []
    skip: set[i64] = set()
    cdef double absolute_rate = top.absolute_rate or -1.

    for i, frag_i in top.frag_list.items():
        if i in skip:
            continue
        uni_rx = top.reactions.get((frag_i.name, None, None, None))
        if uni_rx is not None and len(uni_rx) > 0:
            for rx in uni_rx:
                if detection_one([frag_i], rx, pos, box, absolute_rate):
                    skip.add(i)
                    reactions.append(([frag_i], rx))
                    break
            if i in skip:
                continue
        if not top.r_continue.get((frag_i.name, None, None)):
            continue
        atomi, _ = top.neighbor_atom_map[frag_i.name]
        part_id = frag_i.atoms[atomi]
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
                    if detection_one([frag_i, frag_j], rx, pos, box, absolute_rate):
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
                        if detection_one([frag_i, frag_j, frag_k], rx, pos, box, absolute_rate):
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
                            if detection_one([frag_i, frag_j, frag_k, frag_l], rx, pos, box, absolute_rate):
                                skip.add(i)
                                skip.add(j)
                                skip.add(k)
                                skip.add(l)
                                reactions.append(([frag_i, frag_j, frag_k, frag_l], rx))
                                break

    return reactions

