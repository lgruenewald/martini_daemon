use numpy::{PyReadwriteArray2, PyUntypedArrayMethods};
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use std::collections::{HashSet};
use pyo3_stub_gen::{derive::gen_stub_pyclass, derive::gen_stub_pymethods};


use crate::periodic_box::PeriodicBox;

#[gen_stub_pyclass]
#[pyclass]
pub struct BondGraph {
    n_atoms: usize,
    bonds: Vec<HashSet<usize>>
}

#[gen_stub_pymethods]
#[pymethods]
impl BondGraph {
    #[new]
    /// Create a new empty bond graph.
    ///
    /// Currently bond graphs are internally Vec<HashSet<usize>>
    pub fn new(n_atoms: usize) -> Self {
        BondGraph {
            n_atoms,
            bonds: (0..n_atoms).map(|_| HashSet::new()).collect()
        }
    }

    /// Add a bond to the bond graph.
    ///
    /// Note: i==j, or already present bonds will be ignored.
    ///
    /// i<j or j>i does not matter.
    pub fn add_bond(&mut self, i: usize, j: usize) {
        if i != j {
            self.bonds[i].insert(j);
            self.bonds[j].insert(i);
        }
    }

    /// Make positions whole (in place).
    ///
    /// Uses bond graph in self, and pbc passed as first argument.
    ///
    /// Performs a depth-first traversal of the bond graph and mutates pos, so each bond discovered
    /// is made as short as possible by translating by whole multiples of the periodic box vectors.
    pub fn make_whole<'py>(&self, pbc: &PeriodicBox, mut pos: PyReadwriteArray2<f64>) -> PyResult<()> {
        if pos.shape()[0] != self.n_atoms {
            return Err(PyValueError::new_err("Supplied positions have wrong dimension. Is the number of atoms correct?"));
        }
        if pos.shape()[1] != 3 {
            return Err(PyValueError::new_err("Supplied positions second dimension is not 3."))
        }
        let pos = pos.as_slice_mut()?;

        let mut visited: Vec<bool> = (0..self.n_atoms).map(|_| false).collect();
        for i in 0..self.n_atoms {
            if visited[i] {
                continue;
            }

            let mut stack: Vec<(usize, [f64; 3])> = Vec::new();
            stack.push((i, [pos[i*3 + 0], pos[i*3 + 1], pos[i*3 + 2]]));

            while stack.len() > 0 {
                let (c, reference) = stack.pop().unwrap();
                if visited[c] {
                    continue;
                }
                let cpos: [f64; 3] = pos[c*3..(c+1)*3].try_into()?;
                pos[c*3..(c+1)*3].copy_from_slice(&pbc.move_to(reference, cpos));

                visited[c] = true;

                for j in self.bonds[c].iter() {
                    stack.push((*j, pos[c*3..(c+1)*3].try_into()?));
                }
            }
        }
        for v in visited {
            assert!(v);
        }

        Ok(())
    }

    /// Get a list of bonds.
    ///
    /// Each combination (i, j) is guaranteed to only be there once.
    pub fn to_list(&self) -> Vec<(usize, usize)> {
        let mut res: Vec<(usize, usize)> = Vec::new();
        for i in 0..self.n_atoms {
            for j in self.bonds[i].iter() {
                if i < *j {
                    res.push((i, *j));
                }
            }
        }

        res
    }

    /// Return which atoms can be reached by graph traversal.
    ///
    /// Will perform a depth-first traversal, and include the starting atoms in the result.
    pub fn reachable_from(&self, atoms: HashSet<usize>) -> HashSet<usize> {
        let mut res: HashSet<usize> = HashSet::new();
        let mut stack: Vec<usize> = atoms.into_iter().collect();
        while stack.len() > 0 {
            let atom = stack.pop().unwrap();
            if res.contains(&atom) {
                continue;
            }
            res.insert(atom);
            for j in self.bonds[atom].iter() {
                stack.push(*j);
            }
        }
        res
    }
}
