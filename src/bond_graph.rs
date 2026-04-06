use numpy::{PyArray2, PyArrayMethods, PyUntypedArrayMethods};
use petgraph::data::FromElements;
use pyo3::exceptions::PyValueError;
use pyo3::{prelude::*};
use petgraph::prelude::*;
use petgraph::algo::min_spanning_tree;

use crate::periodic_box::PeriodicBox;

#[pyclass]
pub struct BondGraph {
    n_atoms: usize,
    bonds: UnGraph<(), ()>
}

#[pymethods]
impl BondGraph {
    #[new]
    pub fn new(n_atoms: usize) -> Self {
        let mut bonds = UnGraph::new_undirected();
        for _ in 0..n_atoms {
            bonds.add_node(());
        }
        BondGraph {
            n_atoms,
            bonds
        }
    }

    pub fn add_bond(&mut self, i: u32, j: u32) {
        self.bonds.extend_with_edges(&[(i, j)]);
    }

    pub fn make_whole<'py>(&self, pbc: &PeriodicBox, pos: Bound<'py, PyArray2<f64>>) -> PyResult<()> {
        if pos.shape()[0] != self.n_atoms {
            return Err(PyValueError::new_err("Supplied positions have wrong dimension. Is the number of atoms correct?"));
        }
        if pos.shape()[1] != 3 {
            return Err(PyValueError::new_err("Supplied positions second dimension is not 3."))
        }

        // 1. build a minimum spanning tree from the bond graph
        let mut arr = unsafe { pos.as_array_mut() };
        let mst = UnGraph::from_elements(min_spanning_tree(&self.bonds));

        // 2. use move_to to make all trees whole across PBC
        let mut visited: Vec<bool> = (0..self.n_atoms).map(|_| false).collect();
        for i in 0..self.n_atoms {
            if visited[i] {
                continue;
            }
            visited[i] = true;
            let reference: [f64; 3] = [arr[[i, 0]], arr[[i, 1]], arr[[i, 2]]];
            let mut dfs = Dfs::new(&mst, NodeIndex::<u32>::new(i));
            while let Some(j) = dfs.next(&mst) {
                let j = j.index();
                if i == j {
                    continue;
                }
                assert!(!visited[j]);
                let before = [arr[[j, 0]], arr[[j, 1]], arr[[j, 2]]];
                let after = pbc.move_to(reference, before);
                arr[[j, 0]] = after[0];
                arr[[j, 1]] = after[1];
                arr[[j, 2]] = after[2];
                visited[j] = true;
            }
        }
        for v in visited {
            assert!(v);
        }

        Ok(())
    }

    pub fn to_list(&self) -> Vec<(usize, usize)> {
        let mut res: Vec<(usize, usize)> = Vec::new();

        for edge in self.bonds.raw_edges() {
            let i = edge.source().index();
            let j = edge.target().index();
            res.push((i, j));
        }

        res
    }
}
