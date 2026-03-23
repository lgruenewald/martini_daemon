use std::collections::HashSet;
use numpy::{Ix2, PyArray};
use pyo3::prelude::*;
use rand::prelude::*;

use crate::{
    fragment::Fragment,
    detection_template::DetectionTemplate,
    periodic_box::PeriodicBox,
};

pub fn detection_one<'py>(
    mut rx: PyRefMut<'py, DetectionTemplate>,
    frags: Vec<&Fragment>,
    pbc: &PeriodicBox,
    pos: Borrowed<PyArray<f64, Ix2>>,
    absolute_rate: Option<f64>,
    rng: &mut ThreadRng,
) -> bool {

    // check max distances
    for (idi, atom_i, idj, atom_j, max) in rx.distance_max.iter() {
        let p1 = frags[*idi].atoms[*atom_i];
        let p2 = frags[*idj].atoms[*atom_j];
        if p1 == -1 || p2 == -1 {
            // missing optional atoms
            continue;
        }
        let pos1 = pos
            .get_item(p1)
            .expect("p1 out of range")
            .extract::<[f64; 3]>()
            .expect("can't extract [f64; 3]");
        let pos2 = pos
            .get_item(p2)
            .expect("p2 out of range")
            .extract::<[f64; 3]>()
            .expect("can't extract [f64; 3]");
        let dist = pbc.distance_squared(pos1, pos2);
        if dist > (*max * *max) {
            return false;
        }
    }

    // check min distances
    for (idi, atom_i, idj, atom_j, min) in rx.distance_min.iter() {
        let p1 = frags[*idi].atoms[*atom_i];
        let p2 = frags[*idj].atoms[*atom_j];
        if p1 == -1 || p2 == -1 {
            // missing optional atoms
            continue;
        }
        let pos1 = pos
            .get_item(p1)
            .expect("p1 out of range")
            .extract::<[f64; 3]>()
            .expect("can't extract [f64; 3]");
        let pos2 = pos
            .get_item(p2)
            .expect("p2 out of range")
            .extract::<[f64; 3]>()
            .expect("can't extract [f64; 3]");
        let dist = pbc.distance_squared(pos1, pos2);
        if dist < (*min * *min) {
            return false;
        }
    }

    // check for overlapping atoms
    let mut set: HashSet<isize> = HashSet::new();
    for frag in frags.iter() {
        for p in frag.atoms.iter() {
            if *p == -1 {
                continue;
            }
            if set.contains(p) {
                return false;
            }
            set.insert(*p);
        }
    }

    // check angles
    for (idi, atom_i, idj, atom_j, idk, atom_k, min, max) in rx.angle_limits.iter() {
        let p1 = frags[*idi].atoms[*atom_i];
        let p2 = frags[*idj].atoms[*atom_j];
        let p3 = frags[*idk].atoms[*atom_k];

        if p1 == -1 || p2 == -1 || p3 == -1 {
            // missing optional atoms
            continue;
        }
        let pos1 = pos
            .get_item(p1)
            .expect("p1 out of range")
            .extract::<[f64; 3]>()
            .expect("can't extract [f64; 3]");
        let pos2 = pos
            .get_item(p2)
            .expect("p2 out of range")
            .extract::<[f64; 3]>()
            .expect("can't extract [f64; 3]");
        let pos3 = pos
            .get_item(p3)
            .expect("p3 out of range")
            .extract::<[f64; 3]>()
            .expect("can't extract [f64; 3]");
        let cos_angle = pbc.cos_angle(pos1, pos2, pos3);
        if cos_angle <= *min && cos_angle >= *max {
            // inverted comparison
            // cos is constantly decreasing
            // cos_min => min angle => max cosine value
            // cos_max => max angle => min cosine value
            return false;
        }
    }

    // check dihedrals
    for (idi, atom_i, idj, atom_j, idk, atom_k, idl, atom_l, min, max) in rx.dihedral_limits.iter() {
        let p1 = frags[*idi].atoms[*atom_i];
        let p2 = frags[*idj].atoms[*atom_j];
        let p3 = frags[*idk].atoms[*atom_k];
        let p4 = frags[*idl].atoms[*atom_l];

        if p1 == -1 || p2 == -1 || p3 == -1 || p4 == -1 {
            // missing optional atoms
            continue;
        }
        let pos1 = pos
            .get_item(p1)
            .expect("p1 out of range")
            .extract::<[f64; 3]>()
            .expect("can't extract [f64; 3]");
        let pos2 = pos
            .get_item(p2)
            .expect("p2 out of range")
            .extract::<[f64; 3]>()
            .expect("can't extract [f64; 3]");
        let pos3 = pos
            .get_item(p3)
            .expect("p3 out of range")
            .extract::<[f64; 3]>()
            .expect("can't extract [f64; 3]");
        let pos4 = pos
            .get_item(p4)
            .expect("p4 out of range")
            .extract::<[f64; 3]>()
            .expect("can't extract [f64; 3]");
        let angle = pbc.dihedral(pos1, pos2, pos3, pos4);
        if angle >= *min && angle <= *max {
            return false;
        }
    }

    // simple probability based rate control
    assert!(rx.probability <= 1.);
    if rx.probability < 1. {
        assert!(rx.probability >= 0.);
        let rand: f64 = rng.random();
        assert!(rand >= 0.);
        assert!(rand <= 1.);
        if rand > rx.probability {
            return false;
        }
    }

    // is reaction rate controlled
    if let Some(_) = rx.relative_rate {
        rx.reaction_counter += 1;
        let Some(abs_rate) = absolute_rate else {
            // rate controlled reaction but no global rate control established yet
            return false;
        };
        let Some(obs_rate) = rx.observed_rate else {
            // warmup for this reaction
            return false;
        };
        let prob = abs_rate / obs_rate;
        assert!(prob >= 0.);
        assert!(prob <= 1.);
        let rand: f64 = rng.random();
        assert!(rand >= 0.);
        assert!(rand <= 1.);
        if prob < rand {
            // random between 0 and 1
            return false;
        }
    }

    true
}