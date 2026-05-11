use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3_stub_gen::{derive::gen_stub_pyclass, derive::gen_stub_pymethods};
use std::collections::HashSet;

#[gen_stub_pyclass]
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
    pub angle_limits: Vec<(usize, usize, usize, usize, usize, usize, f64, f64)>,
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
    #[pyo3(get, set)]
    pub probability: f64,
}

#[gen_stub_pymethods]
#[pymethods]
impl DetectionTemplate {
    #[new]
    /// Create a new empty detection template.
    ///
    /// Must be fully parsed before adding it to DetectionTemplateList.
    pub fn new() -> Self {
        Self {
            name: None,
            reactants: Vec::new(),
            distance_max: Vec::new(),
            distance_min: Vec::new(),
            angle_limits: Vec::new(),
            dihedral_limits: Vec::new(),
            probability: 1.,
        }
    }

    /// Add a maximum distance condition.
    ///
    /// Entry should be: reactant idx, graph node idx, reactant idx, graph node idx, distance^2 (nm^2)
    pub fn add_distance_max(&mut self, entry: (usize, usize, usize, usize, f64)) {
        self.distance_max.push(entry);
    }

    /// Add a minimum distance condition.
    ///
    /// Entry should be: reactant idx, graph node idx, reactant idx, graph node idx, distance^2 (nm^2)
    pub fn add_distance_min(&mut self, entry: (usize, usize, usize, usize, f64)) {
        self.distance_min.push(entry);
    }

    /// Add a forbidden angle range.
    ///
    /// Entry should be: reactant idx, graph node (3x), minimum (cosine of angle), maximum (cosine of angle)
    pub fn add_angle_limit(&mut self, entry: (usize, usize, usize, usize, usize, usize, f64, f64)) {
        self.angle_limits.push(entry);
    }

    /// Add a forbidden dihedral range.
    ///
    /// Entry should be: reactant idx, graph node (4x), minimum (angle), maximum (angle).
    ///
    /// min < max must be true. both min and max must be between -pi and pi, in radians.
    pub fn add_dihedral_limit(
        &mut self,
        #[gen_stub(override_type(
            type_repr = "tuple[int, int, int, int, int, int, int, int, float, float]"
        ))]
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

    /// Get a list of angle limits already present.
    ///
    /// Used during parsing.
    #[getter]
    #[gen_stub(override_return_type(
        type_repr = "list[tuple[int, int, int, int, int, int, float, float]]"
    ))]
    pub fn angle_limits(&self) -> Vec<(usize, usize, usize, usize, usize, usize, f64, f64)> {
        self.angle_limits.clone()
    }

    /// Get a list of dihedral limits already present.
    ///
    /// Used during parsing.
    #[getter]
    #[gen_stub(override_return_type(
        type_repr = "list[tuple[int, int, int, int, int, int, int, int, float, float]]"
    ))]
    pub fn dihedral_limits(
        &self,
    ) -> Vec<(
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
    )> {
        self.dihedral_limits.clone()
    }

    /// Raise an exception if reaction is not valid.
    ///
    /// Must be called when reaction is done parsing.
    pub fn complete(&self) -> Result<(), PyErr> {
        // 1. must have name and reactants
        if self.name.is_none() {
            return Err(PyException::new_err("Reaction must have a name"));
        }
        if self.reactants.is_empty() {
            return Err(PyException::new_err(format!(
                "Reaction {} has no reactants.",
                self.name.clone().unwrap()
            )));
        }
        // 2. must be connected through r_max
        let mut connections: Vec<HashSet<usize>> =
            self.reactants.iter().map(|_| HashSet::new()).collect();
        for (m1, _i1, m2, _i2, _dist) in self.distance_max.iter() {
            // parse_pair() in the parsing code should already handle out of bounds
            connections[*m1].insert(*m2);
            connections[*m2].insert(*m1);
        }

        let mut marked: HashSet<usize> = HashSet::new();
        let mut stack: Vec<usize> = Vec::new();
        stack.push(0);

        while let Some(cur) = stack.pop() {
            marked.insert(cur);
            for i in connections[cur].iter() {
                if !marked.contains(i) {
                    stack.push(*i);
                }
            }
        }

        let all: HashSet<usize> = (0..self.reactants.len()).collect();
        if let Some(unconnected) = marked.difference(&all).next() {
            return Err(PyException::new_err(format!(
                "Reaction {}: Reactants not all connected by r_max. First disconnected reactant index {}.",
                self.name.clone().unwrap(),
                *unconnected
            )));
        }
        Ok(())
    }
}
