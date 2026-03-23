use std::collections::HashSet;
use ordermap::OrderSet;
use numpy::{Ix2, PyArray};
use pyo3::prelude::*;
use kdtree::{KdTree, distance::squared_euclidean};
use sortedlist_rs::SortedList;

use crate::detection_template::DetectionTemplate;
use crate::detection_template_list::DetectionTemplateList;
use crate::fragment::Fragment;
use crate::detection_one::detection_one;
use crate::periodic_box::PeriodicBox;
use crate::frag_list::FragList;

/// Main entry point for the detection algorithm
/// - takes T* components FragList and DetectionTemplateList, a pbc and the current positions.
/// - builds a current frame neighbor list.
/// - runs the detection algorithm on all possible reaction-reactant combinations.
/// - returns a list of reactions, each as a list of frag_ids and reaction names.
#[pyfunction]
pub fn detection<'py>(
    frag_list: &mut FragList,
    detection_template_list: &mut DetectionTemplateList,
    absolute_rate: Option<f64>,
    pbc: &PeriodicBox,
    pos: Bound<'py, PyArray<f64, Ix2>>,
) -> PyResult<Vec<(Vec<usize>, String)>> {
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
        for frag_atom_index in detection_template_list.frag_name_to_r_max_atoms[name].iter() {
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
    let mut skip: OrderSet<usize> = OrderSet::new();
    // result
    let mut reactions: Vec<(Vec<usize>, String)> = Vec::new();

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
                let rx = detection_template_list.reactions.get(uni_rx).unwrap().borrow_mut(py);;
                let frags = vec![frag_i];
                if detection_one(rx, frags, pbc, pos, absolute_rate, &mut rng) {
                    skip.insert(*i);
                    reactions.push((vec![*i], uni_rx.clone()));
                    break; // uni_rx in uni_rxs
                }
            }
        }

        // reacted => skip
        if skip.contains(i) {
            continue;
        }

        // no bimolecular? => skip
        if detection_template_list.first_reactants.get(&frag_i.name).is_none() {
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
                matches.insert(*new_match);
            }
        }

        // TODO tri, tetramolecular reactions
        let mut next_j = 0;
        for j in matches.iter() {
            if *j < next_j {
                continue;
            }
            next_j = *j + 1;

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
                        .borrow_mut(py);
                    let frags = vec![frag_i, frag_j];

                    if detection_one(rx, frags, pbc, pos, absolute_rate, &mut rng) {
                        skip.insert(*i);
                        skip.insert(*j);
                        reactions.push((vec![*i, *j], bi_rx.clone()));
                        break; // bi_rx in bi_rxs
                    }
                }
            };
        }
    }

    Ok(reactions)
}
