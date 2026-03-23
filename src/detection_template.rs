use pyo3::{prelude::*};
use std::{collections::HashSet};
use pyo3::exceptions::PyException;

#[pyclass]
pub struct DetectionTemplate {
    #[pyo3(get, set)]
    pub name: Option<String>,
    #[pyo3(get, set)]
    pub reactants: Vec<String>,
    #[pyo3(get)]
    pub distance_max: Vec<(usize, usize, usize, usize, f64)>,
    #[pyo3(get)]
    pub distance_min: Vec<(usize, usize, usize, usize, f64)>,
    #[pyo3(get)]
    pub angle_limits: Vec<(usize, usize, usize, usize, usize, usize, f64, f64)>,
    #[pyo3(get)]
    pub dihedral_limits: Vec<(
        usize,
        usize,
        usize,
        usize,
        usize,
        usize,
        usize,
        usize,
        f64,
        f64,
    )>,
    // rate stuff
    #[pyo3(get, set)]
    pub relative_rate: Option<f64>,
    #[pyo3(get, set)]
    pub reaction_counter: usize,
    #[pyo3(get, set)]
    pub observed_rate: Option<f64>,
    #[pyo3(get, set)]
    pub probability: f64,
}

#[pymethods]
impl DetectionTemplate {
    #[new]
    pub fn new() -> Self {
        Self {
            name: None,
            reactants: Vec::new(),
            distance_max: Vec::new(),
            distance_min: Vec::new(),
            angle_limits: Vec::new(),
            dihedral_limits: Vec::new(),
            relative_rate: None,
            reaction_counter: 0,
            observed_rate: None,
            probability: 1.,
        }
    }

    pub fn add_distance_max(
        &mut self,
        entry: (usize, usize, usize, usize, f64),
    ) {
        self.distance_max.push(entry);
    }

    pub fn add_distance_min(
        &mut self,
        entry: (usize, usize, usize, usize, f64),
    ) {
        self.distance_min.push(entry);
    }

    pub fn add_angle_limit(
        &mut self,
        entry: (usize, usize, usize, usize, usize, usize, f64, f64),
    ) {
        self.angle_limits.push(entry);
    }

    pub fn add_dihedral_limit<'py>(
        &mut self,
        entry: (
            usize,
            usize,
            usize,
            usize,
            usize,
            usize,
            usize,
            usize,
            f64,
            f64,
        ),
    ) {
        self.dihedral_limits.push(entry);
    }

    /// Raises an exception if reaction is not valid.
    /// Must be called when reaction is done parsing.
    pub fn complete(&self) -> Result<(), PyErr> {
        // 1. must have name and reactants
        if self.name.is_none() {
            return Err(PyException::new_err("Reaction must have a name"))
        }
        if self.reactants.len() == 0 {
            return Err(PyException::new_err(format!(
                "Reaction {} has no reactants.",
                self.name.clone().unwrap()
            )));
        }
        // 2. must be connected through r_max
        let mut connections: Vec<HashSet<usize>> =
            self.reactants.iter().map(|_| HashSet::new()).collect();
        for (m1, i1, m2, i2, dist) in self.distance_max.iter() {
            // parse_pair() in the parsing code should already handle out of bounds
            connections[*m1].insert(*m2);
            connections[*m2].insert(*m1);
        }

        let mut marked: HashSet<usize> = HashSet::new();
        let mut stack: Vec<usize> = Vec::new();
        stack.push(0);

        while stack.len() > 0 {
            let cur = stack.pop().unwrap();
            marked.insert(cur);
            for i in connections[cur].iter() {
                if !marked.contains(i) {
                    stack.push(*i);
                }
            }
        }

        let all: HashSet<usize> = (0..self.reactants.len()).collect();
        for unconnected in marked.difference(&all) {
            return Err(PyException::new_err(format!(
                "Reaction {}: Reactants not all connected by r_max. First disconnected reactant index {}.",
                self.name.clone().unwrap(), *unconnected
            )));
        }
        Ok(())
    }
}
