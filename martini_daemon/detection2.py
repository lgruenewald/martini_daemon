# detection algo, attempt 2
# NOTE: currently broken
import numpy as np
from numpy.random import Generator
from .reaction_template import ReactionTemplate
from .fragment import Fragment
from freud.box import Box
from freud.locality import AABBQuery


def pbc_diff(p1, p2, size):
    """
    Difference (p1-p2) of two positions in a pbc box of size.
    """
    diff = p1 - p2
    base = np.floor(diff / size + 0.5) * size
    return diff - base


def pbc_dist(p1, p2, size):
    """
    Distance between two points in a pbc box of size.
    """
    return np.sqrt(np.sum(pbc_diff(p1, p2, size)**2, axis=-1))


def pbc_angle_cos(p1, p2, p3, size):
    """
    Cosine of the angle enclosed by p1, p2 and p3, p2 being the central atom,
    in a pbc box of size.
    """
    v1 = pbc_diff(p2, p1, size)
    v2 = pbc_diff(p2, p3, size)
    return (
        np.vecdot(v1, v2)
        / np.linalg.norm(v1, axis=-1)
        / np.linalg.norm(v2, axis=-1)
    )


def pbc_dihedral(p1, p2, p3, p4, size):
    """
    Dihedral (as defined by IUPAC) between the atoms at p1, p2, p3, p4,
    in a pbc box of size. Returns it between 0 and 2 pi (in radians).
    """
    # difference vectors
    u1 = pbc_diff(p2, p1, size)
    u2 = pbc_diff(p3, p2, size)
    u3 = pbc_diff(p4, p3, size)
    # cross products
    v1 = np.cross(u1, u2)
    v2 = np.cross(u2, u3)
    v3 = np.cross(v1, v2)
    # x and y
    x = np.vecdot(v1, v2)
    y = np.vecdot(u2, v3) / np.linalg.norm(u2, axis=-1)
    return np.mod(
        np.atan2(y, x),
        np.pi * 2
    )


def pbc_crosses(p1, p2, size):
    return np.any(
        np.abs(p1 - p2) / size > 0.5,
        axis=-1
    )


def detection_all(
    all_fragments: dict[int, Fragment],
    all_templates: list[ReactionTemplate],
    box, pos: np.ndarray,
    nlist_cutoff: float,
    absolute_rate: float,
    generator: Generator
) -> list[tuple[list[Fragment], ReactionTemplate]]:

    reactions: list[list[Fragment], ReactionTemplate] = []

    if len(all_fragments) == 0 or len(all_templates) == 0:
        return reactions

    # preparation for detection
    # 1. make lists of fragments of each type
    # TODO if this works change the T* storage of fragments to make
    # this construction here easier, possibly already numpy-ified?
    # relying on the fact that python dicts preserve insertion order
    frag_types: dict[str, list[Fragment]] = {}
    frag_names: dict[str, int] = {}  # frag name -> frag name index
    # split all fragments into a per type list of fragments
    for frag_id, frag in all_fragments.items():
        if frag_types.get(frag.name) is None:
            frag_names[frag.name] = len(frag_names)
            frag_types[frag.name] = []
        frag_types[frag.name].append(frag)
    # make a numpy representation of the fragment list
    frag_data: list[np.ndarray] = []
    # where each fragment begins/ends, if frag_data was contiguous
    # this is how they are ordered for the AABBQuery
    offsets: list[int] = [0]
    for name, frags in frag_types.items():
        atoms_per_frag = len(frags[0].atoms)
        # 2D array for each frag type
        # columns: frag_id atom1 atom2 ... atomn
        data = np.empty(
            (len(frags), atoms_per_frag + 1),
            dtype=np.int64
        )
        data[:, 0] = [frag.frag_id for frag in frags]
        data[:, 1:] = [frag.atoms for frag in frags]
        frag_data.append(data)
        offsets.append(offsets[-1] + len(frags))

    # 2. build Axis aligned bounding box (AABB) query for all fragments once
    pos_filtered = np.empty((len(all_fragments), 3), pos.dtype)
    for i, (name, i2) in enumerate(frag_names.items()):
        assert i == i2
        start = offsets[i]
        end = offsets[i+1]
        # currently, 1st atom must be a normal atom (see graph.py validation)
        pos_filtered[start:end] = pos[frag_data[i][:, 1]]
    query = AABBQuery(
        Box(box[0], box[1], box[2]),
        pos_filtered
    )

    # 4. skipset of frag_ids
    # yeah... this is the main obstacle for good parallelization
    skipset = set()

    for rx in all_templates:
        if len(rx.reactants) > 2:
            # TODO temporary limitation for performance
            raise ValueError("More than 2 reactants currently not supported.")
        # for each first reactant:
        # fr short for First Reactant
        fr = rx.reactants[0]
        if frag_names.get(fr) is None:
            continue
        fr_i = frag_names[fr]
        fr_start = offsets[fr_i]
        fr_count = len(frag_data[fr_i])

        fr_frag_ids = frag_data[fr_i][:, 0]
        fr_atoms = frag_data[fr_i][:, 1:]
        fr_pos = pos[frag_data[fr_i][:, 1:]]
        # lr short for Last Reactant
        if len(rx.reactants) > 1:
            lr = rx.reactants[1]
            if frag_names.get(lr) is None:
                continue
            lr_i = frag_names[lr]
            lr_start = offsets[lr_i]
            lr_end = offsets[lr_i + 1]

        for i in range(fr_count):
            # 1. continue if first reactant is in the skiplist
            fr_frag_id = fr_frag_ids[i]
            fr_frag_i = fr_start + i
            if fr_frag_id in skipset:
                continue

            # 2. take the first reactant, build a neighbor list around it,
            # or use the cached neighbor list around it if it was already built
            # TODO: when more than 2 reactants, make a neighbor list union
            # sort of like tuples of combinations within a distance:
            # (1, 2), (1, 3), (2, 3) + a dense array of candidates for the last
            # reactant and do the geometry check vectorized on the list for
            # the last reactant. Hopefully there is a way to make the "union"
            # fast? because performance will also depend on that
            if len(rx.reactants) > 1:
                nlist = query.query(
                    fr_pos[i],
                    {"r_max": nlist_cutoff}
                ).toNeighborList()
                nlist.filter(
                    (nlist[:, 1] >= lr_start) & (nlist[:, 1] < lr_end)
                    # no bimolecular self reaction
                    & (nlist[:, 1] != fr_frag_i)
                )
                lr_list = nlist[:, 1] - lr_start

                # 3. build a numpy array with coordinates for neighbors,
                # of frag_id (shape: (lr_list_len))
                lr_frag_ids = frag_data[lr_i][lr_list, 0]
                # of list of atoms (shape: (lr_list_len, atoms_per_frag))
                lr_atoms = frag_data[lr_i][lr_list, 1:]
                # and of positions (shape: (lr_list_len, atoms_per_frag, 3))
                lr_pos = pos[lr_atoms]
                # which in lr pass
                reaction_pass: np.ndarray = np.repeat(True, len(lr_frag_ids))
            else:
                reaction_pass = np.array([True])

            # 4. run the geom conditions vectorized
            # - further complexity - optional atoms can be missing

            if len(rx.reactants) == 1:
                atoms = [fr_atoms[i:i+1]]
                positions = [fr_pos[i:i+1]]
            else:
                atoms = [fr_atoms[i:i+1], lr_atoms]
                positions = [fr_pos[i:i+1], lr_pos]

            #  distance checks
            for idi, atomi, idj, atomj, max in rx.distance_max:
                # missing optional atom
                missing_opt = (
                    (atoms[idi][:, atomi] == -1)
                    | (atoms[idj][:, atomj] == -1)
                )
                within_cutoff = pbc_dist(
                    positions[idi][:, atomi, :],
                    positions[idj][:, atomj, :],
                    box
                ) < max
                reaction_pass = reaction_pass | (missing_opt & within_cutoff)
            for idi, atomi, idj, atomj, min in rx.distance_min:
                # missing optional atom
                missing_opt = (
                    (atoms[idi][:, atomi] == -1)
                    | (atoms[idj][:, atomj] == -1)
                )
                within_cutoff = pbc_dist(
                    positions[idi][:, atomi, :],
                    positions[idj][:, atomj, :],
                    box
                ) > min
                reaction_pass = reaction_pass | (missing_opt & within_cutoff)

            for (
                idi, atomi, idj, atomj, idk, atomk, min, max
            ) in rx.angle_limits:
                missing_opt = (
                    (atoms[idi][:, atomi] == -1)
                    | (atoms[idj][:, atomj] == -1)
                    | (atoms[idk][:, atomk] == -1)
                )
                angle_cos = pbc_angle_cos(
                    positions[idi][:, atomi, :],
                    positions[idj][:, atomj, :],
                    positions[idk][:, atomk, :],
                    box
                )
                within_cutoff = (angle_cos > min) & (angle_cos < max)
                reaction_pass = reaction_pass | (missing_opt & within_cutoff)

            for (
                idi, atomi, idj, atomj, idk, atomk, idl, atoml, min, max
            ) in rx.dihedral_limits:
                missing_opt = (
                    (atoms[idi][:, atomi] == -1)
                    | (atoms[idj][:, atomj] == -1)
                    | (atoms[idk][:, atomk] == -1)
                    | (atoms[idl][:, atoml] == -1)
                )
                dihedral = pbc_dihedral(
                    positions[idi][:, atomi, :],
                    positions[idj][:, atomj, :],
                    positions[idk][:, atomk, :],
                    positions[idl][:, atoml, :],
                    box
                )
                within_cutoff = (dihedral > min) & (dihedral < max)
                reaction_pass = reaction_pass | (missing_opt & within_cutoff)

            # 5. run the rate condition vectorized
            if rx.relative_rate is not None:
                rx.reaction_counter = reaction_pass.sum()
                if absolute_rate <= 0. or rx.observed_rate <= 0.:
                    # negative/0 absolute or observed rate means warmup
                    reaction_pass[:] = False
                else:
                    # prob of 1.0 = always, 0.0 = never
                    prob = absolute_rate / rx.observed_rate
                    reaction_pass = reaction_pass & (
                        generator.random(reaction_pass.shape)
                        < prob
                    )

            # 6. filter reactions that can happen
            # we only need to find the first reaction that fulfills
            # these conditions at this point
            if len(rx.reactants) > 1:
                candidates = []
                for j in range(len(lr_frag_ids)):
                    lr_frag_id = lr_frag_ids[j]
                    # - must have passed geom+rate conditions
                    if not reaction_pass[j]:
                        continue
                    # - every fragment can take part only in 1 reaction
                    if lr_frag_id in skipset:
                        continue
                    # - overlapping fragments can't react
                    overlap = np.intersect1d(
                        fr_atoms,
                        lr_atoms[j]
                    )
                    if (
                        len(overlap) == 0
                        or len(overlap) == 1 and overlap[0] == -1
                    ):
                        continue
                    candidates.append(lr_frag_id)

                # 7. modify state
                if len(candidates) > 0:
                    candidate = generator.choice(candidates)
                    # - update the skipset in this process
                    skipset.add(fr_frag_id)
                    skipset.add(candidate)

                    # - add valid reaction to reactions
                    reactions.append((
                        [
                            all_fragments[fr_frag_id],
                            all_fragments[candidate]
                        ],
                        rx
                    ))
            else:
                # unimolecular
                if reaction_pass[0]:
                    skipset.add(fr_frag_id)
                    reactions.append((
                        [all_fragments[fr_frag_id]],
                        rx
                    ))

    return reactions
