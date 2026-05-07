use std::collections::HashSet;
use numpy::{Ix2, PyArray};
use pyo3::prelude::*;
use kdtree::{KdTree, distance::squared_euclidean};
use sortedlist_rs::SortedList;
use pyo3_stub_gen::{derive::gen_stub_pyfunction};

use crate::detection_template_list::DetectionTemplateList;
use crate::detection_one::detection_one;
use crate::periodic_box::PeriodicBox;
use crate::frag_list::FragList;

/// Perform the detection algorithm.
///
/// Main entry point for the detection algorithm
/// - takes T* components FragList and DetectionTemplateList, a pbc and the current positions.
/// - builds a current frame neighbor list.
/// - runs the detection algorithm on all possible reaction-reactant combinations.
/// - returns a list of reactions, each as a list of frag_ids and reaction names.
#[gen_stub_pyfunction]
#[pyfunction]
pub fn detection<'py>(
    #[gen_stub(override_type(type_repr="FragList"))]
    frag_list: &mut FragList,
    #[gen_stub(override_type(type_repr="DetectionTemplateList"))]
    detection_template_list: &mut DetectionTemplateList,
    pbc: &PeriodicBox,
    pos: Bound<'py, PyArray<f64, Ix2>>,
) -> PyResult<Vec<(String, Vec<usize>)>> {
    let py = pos.py();
    let pos = pos.as_borrowed();

    // 1. KD Tree Build
    let mut tree = KdTree::new(3);
    // some safety margin for floating point inaccuracies
    let cutoff = detection_template_list.largest_r_max + 0.1;

    for frag in frag_list.values() {
        let name = &frag.name;
        let frag_id = frag.frag_id;
        // for each within-frag index
        let Some(name_to_rmax) = detection_template_list.frag_name_to_r_max_atoms.get(name) else {
            continue;
        };
        for frag_atom_index in name_to_rmax.iter() {
            // index globally
            let index = frag.atoms[*frag_atom_index];
            if index < 0 {
                // missing optional atom
                continue;
            }
            let pos = pos.get_item(index as usize)?.extract::<[f64; 3]>()?;
            // should be within the middle pbc copy
            let pos = pbc.move_within(pos);

            // add the copies that are within or almost within (wrt cutoff) PBC to KDTree
            // this is because the KDTree library is pbc unaware
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        // Note: dx=0, dy=0, dz=0 is handled here too
                        let new_pos = pbc.translate_by(pos, dx, dy, dz);
                        if pbc.is_almost_inside(new_pos, cutoff) {
                            tree.add(new_pos, frag_id).unwrap();
                        }
                    }
                }
            }
        }
    }

    // 2. Run detection algorithm on all combinations of reactants and reactions
    // initialize rng only once, will get passed to detection_one
    let mut rng = rand::rng();
    // already reacted - skip
    // note: this is only for (theoretical) speedups
    // the modification algorithm already ensures that each fragment can only react in one rx
    // may be removed for parallelization in the future
    // for many systems skip is going to be quite empty or low in # of members
    let mut skip: HashSet<usize> = HashSet::new();
    // result
    let mut reactions: Vec<(String, Vec<usize>)> = Vec::new();

    for (i, frag_i) in frag_list.iter() {
        if skip.contains(i) {
            continue;
        }
        // Unimolecular reactions
        if let Some(uni_rxs) = detection_template_list
            .reactions_by_reactants
            .get(&vec![frag_i.name.to_string()])
        {
            for uni_rx in uni_rxs {
                let rx = detection_template_list.reactions.get(uni_rx).unwrap().borrow(py);
                let frags = vec![frag_i];
                if detection_one(&rx, frags, pbc, pos, &mut rng) {
                    skip.insert(*i);
                    reactions.push((uni_rx.clone(), vec![*i]));
                    break; // uni_rx in uni_rxs
                }
            }
        }

        // reacted => skip
        if skip.contains(i) {
            continue;
        }

        // no multimolecular? => skip
        let Some((bi, tri)) = detection_template_list.first_reactants.get(&frag_i.name) else {
            continue;
        };
        if !*bi && !*tri {
            continue;
        }

        // build a list of close by frag_IDs
        // we need to start from every r_max having atom unfortunately, to guarantee that all
        // matches are still findable
        let mut matches = SortedList::new();
        for frag_atom_index in detection_template_list.frag_name_to_r_max_atoms[&frag_i.name].iter() {
            let index = frag_i.atoms[*frag_atom_index];
            if index == -1 {
                // missing optional atom
                continue;
            }
            let pos: [f64; 3] = pos
                .get_item(index as usize)?
                .extract::<[f64; 3]>()?;
            for (_, new_match) in tree.within_unsorted(&pos, cutoff.powi(2), &squared_euclidean).unwrap() {
                if !matches.contains(new_match) {
                    matches.insert(*new_match);
                }
            }
        }

        // == BIMOLECULAR REACTIONS ==
        if *bi {
            for j in matches.iter() {
                // no self reactions
                if *i == *j {
                    continue;
                }
                // i can become skipped during a multimolecular reaction too
                if skip.contains(j) || skip.contains(i) {
                    continue;
                }
                let frag_j = frag_list.get(*j).unwrap();
                if let Some(bi_rxs) = detection_template_list
                    .reactions_by_reactants
                    .get(&vec![frag_i.name.to_string(), frag_j.name.to_string()])
                {
                    for bi_rx in bi_rxs {
                        let rx = detection_template_list
                            .reactions
                            .get(bi_rx)
                            .unwrap()
                            .borrow(py);
                        let frags = vec![frag_i, frag_j];

                        if detection_one(&rx, frags, pbc, pos, &mut rng) {
                            skip.insert(*i);
                            skip.insert(*j);
                            reactions.push((bi_rx.clone(), vec![*i, *j]));
                            break; // bi_rx in bi_rxs
                        }
                    }
                };
            }
        }
        // == TRIMOLECULAR REACTIONS ==
        // this is some not so nice code...
        // trimolecular reactions will be somewhat slower to check for
        // it is necessary to extend neighbor list matches because of scenarios where
        // frag_i <==r_max==> frag_j <==r_max==> frag_k
        // as well as
        // frag_i <==r_max==> frag_k <==r_max==> frag_j
        // in some of these scenarios either frag_j or frag_k may not be in the neighbor list of frag_i
        // but since detection_template.complete() ensures that all reactants are connected
        // frag_i must be within neighbor list of either frag_j or frag_k
        if *tri {
            // for now this is some code repetition, but I wanted to keep it isolated
            // from bimolecular code, because that one alone is quite ok
            // this can be refactored later...
            for j_or_k in matches.iter() {
                if *j_or_k == *i {
                    continue;
                }
                if skip.contains(i) || skip.contains(j_or_k) {
                    continue;
                }
                let frag_j_or_k = frag_list.get(*j_or_k).unwrap();
                // new extended neighbor list
                let mut tri_matches = SortedList::new();
                for c in matches.iter() {
                    if !tri_matches.contains(c) {
                        tri_matches.insert(*c);
                    }
                }
                for frag_atom_index in detection_template_list.frag_name_to_r_max_atoms[&frag_j_or_k.name].iter() {
                    let index = frag_j_or_k.atoms[*frag_atom_index];
                    if index == -1 {
                        // missing optional atom
                        continue;
                    }
                    let pos: [f64; 3] = pos
                        .get_item(index as usize)?
                        .extract::<[f64; 3]>()?;
                    for (_, new_match) in tree.within_unsorted(&pos, cutoff.powi(2), &squared_euclidean).unwrap() {
                        if !tri_matches.contains(new_match) {
                            tri_matches.insert(*new_match);
                        }
                    }
                }

                for k_or_j in tri_matches.iter() {
                    if *k_or_j == *i || *k_or_j == *j_or_k {
                        continue;
                    }
                    if skip.contains(i) || skip.contains(j_or_k) || skip.contains(k_or_j) {
                        continue;
                    }
                    let frag_k_or_j = frag_list.get(*k_or_j).unwrap();
                    // get the first reaction that's possible
                    if let Some(tri_rxs) = detection_template_list
                        .reactions_by_reactants
                        .get(&vec![frag_i.name.to_string(), frag_j_or_k.name.to_string(), frag_k_or_j.name.to_string()])
                    {
                        for tri_rx in tri_rxs {
                            let rx = detection_template_list
                                .reactions
                                .get(tri_rx)
                                .unwrap()
                                .borrow(py);
                            let frags = vec![frag_i, frag_j_or_k, frag_k_or_j];

                            if detection_one(&rx, frags, pbc, pos, &mut rng) {
                                skip.insert(*i);
                                skip.insert(*j_or_k);
                                skip.insert(*k_or_j);
                                reactions.push((tri_rx.clone(), vec![*i, *j_or_k, *k_or_j]));
                                break; // tri_rx in tri_rxs
                            }
                        }
                    };
                    if skip.contains(i) {
                        break; // for k_or_j
                    }
                    if let Some(tri_rxs) = detection_template_list
                        .reactions_by_reactants
                        .get(&vec![frag_i.name.to_string(), frag_k_or_j.name.to_string(), frag_j_or_k.name.to_string()])
                    {
                        for tri_rx in tri_rxs {
                            let rx = detection_template_list
                                .reactions
                                .get(tri_rx)
                                .unwrap()
                                .borrow(py);
                            let frags = vec![frag_i, frag_k_or_j, frag_j_or_k];

                            if detection_one(&rx, frags, pbc, pos, &mut rng) {
                                skip.insert(*i);
                                skip.insert(*j_or_k);
                                skip.insert(*k_or_j);
                                reactions.push((tri_rx.clone(), vec![*i, *k_or_j, *j_or_k]));
                                break; // tri_rx in tri_rxs
                            }
                        }
                    };


                }



            }
        }
    }

    Ok(reactions)
}
