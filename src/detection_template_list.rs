use ordermap::OrderMap;
use std::collections::HashMap;
use pyo3::prelude::*;
use pyo3::{pyclass, pymethods};
use pyo3::exceptions::PyException;
use crate::detection_template::DetectionTemplate;
#[pyclass]
pub struct DetectionTemplateList {
    pub reactions: OrderMap<String, Py<DetectionTemplate>>,
    pub reactions_by_reactants: OrderMap<Vec<String>, Vec<String>>,
    pub first_reactants: HashMap<String, (bool, bool)>,
    pub frag_name_to_r_max_atoms: HashMap<String, Vec<usize>>,
    pub largest_r_max: f64
}

#[pymethods]
impl DetectionTemplateList {
    #[new]
    pub fn new() -> Self {
        Self {
            reactions: OrderMap::new(),
            reactions_by_reactants: OrderMap::new(),
            first_reactants: HashMap::new(),
            frag_name_to_r_max_atoms: HashMap::new(),
            largest_r_max: 0.0
        }
    }
    pub fn add_detection_template(&mut self, rx: Py<DetectionTemplate>, py: Python) -> PyResult<bool> {
        let rx_ref = rx.borrow(py);

        if rx_ref.reactants.len() > 3 {
            return Err(PyException::new_err("Only up to 3 reactants per reaction are supported."));
        }

        for (i1, a1, i2, a2, r_max) in rx_ref.distance_max.iter() {
            // largest_r_max
            if (i1 == i2) {
                // this is for constructing additional data required for intermolecular
                // reactions and detection algorithm.
                continue;
            }
            // the largest intermolecular r_max, the neighbor list cutoff will be based on this.
            // currently globally a single value, but for most systems this is not a problem
            // as most r_max for simulations will be around ~0.6> nm. The main reason for doing
            // this over hard-coded 1 nm neighbor list cutoff is correctness. Tough,
            // performance should be better in the typical case than the old system with
            // a hardcoded ~1 nm cutoff.
            self.largest_r_max = f64::max(
                self.largest_r_max,
                *r_max
            );

            // frag_name_to_r_max_atoms
            // this stores which atoms have intermolecular r_max for each frag name across
            // all available reactions.
            // this means it's optimal to design all reaction templates for each fragment
            // based on a single r_max on the same atom (this is usually the case!).
            let frag_name1 = &rx_ref.reactants[*i1];
            if !self.frag_name_to_r_max_atoms.contains_key(frag_name1) {
                self.frag_name_to_r_max_atoms.insert(frag_name1.to_string(), Vec::new());
            }
            self.frag_name_to_r_max_atoms.get_mut(frag_name1).unwrap().push(*a1);
            let frag_name2 = &rx_ref.reactants[*i2];
            if !self.frag_name_to_r_max_atoms.contains_key(frag_name2) {
                self.frag_name_to_r_max_atoms.insert(frag_name2.to_string(), Vec::new());
            }
            self.frag_name_to_r_max_atoms.get_mut(frag_name2).unwrap().push(*a2);
        }

        // first reactants - list of reactants that start multimolecular reactions
        if rx_ref.reactants.len() > 1 {
            if !self.first_reactants.contains_key(&rx_ref.reactants[0]) {
                self.first_reactants.insert(rx_ref.reactants[0].clone(), (false, false));
            }
            match (rx_ref.reactants.len(), self.first_reactants.get_mut(&rx_ref.reactants[0])) {
                (2, Some((val, _))) => *val = true,
                (3, Some((_, val))) => *val = true,
                _ => { unreachable!() }
            }

        }
        // reactions by reactants
        let reactants = rx_ref.reactants.clone();
        if let Some(entry) = self.reactions_by_reactants.get_mut(&reactants) {
            entry.push(rx_ref.name.clone().unwrap());
        } else {
            self.reactions_by_reactants.insert(reactants, vec![rx_ref.name.clone().unwrap()]);
        }

        // reactions
        let name = rx_ref.name.clone().unwrap();
        drop(rx_ref);
        self.reactions.insert(name, rx);

        Ok(true)
    }

    pub fn get_detection_template(
        &self,
        name: &str,
        py: Python
    ) -> PyResult<Option<Py<DetectionTemplate>>> {
        // returns a mutable copy of the reaction
        match self.reactions.get(name) {
            Some(r) => Ok(Some(r.clone_ref(py))),
            None => Ok(None),
        }
    }

    pub fn reaction_names(&self) -> Vec<String> {
        self.reactions.keys().map(|a| a.clone()).collect()
    }
}
