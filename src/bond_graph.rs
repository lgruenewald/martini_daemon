use numpy::{PyReadwriteArray2, PyUntypedArrayMethods};
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use std::collections::{HashSet};

use crate::periodic_box::PeriodicBox;

#[pyclass]
pub struct BondGraph {
    n_atoms: usize,
    bonds: Vec<HashSet<usize>>
}

#[pymethods]
impl BondGraph {
    #[new]
    pub fn new(n_atoms: usize) -> Self {
        BondGraph {
            n_atoms,
            bonds: (0..n_atoms).map(|_| HashSet::new()).collect()
        }
    }

    pub fn add_bond(&mut self, i: usize, j: usize) {
        if i != j {
            self.bonds[i].insert(j);
            self.bonds[j].insert(i);
        }
    }

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
