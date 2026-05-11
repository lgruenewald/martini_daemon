use crate::fragment::Fragment;
use ordermap::OrderMap;
use pyo3::prelude::*;
use pyo3_stub_gen::{derive::gen_stub_pyclass, derive::gen_stub_pymethods};
use std::collections::HashMap;
use std::iter::Iterator;
use std::ops::IndexMut;

#[gen_stub_pyclass]
#[pyclass]
pub struct FragList {
    defrag_list: Vec<Vec<usize>>,
    frag_list: OrderMap<usize, Fragment>,
    #[pyo3(get)]
    next_frag_id: usize,
    #[pyo3(get)]
    frag_counts: HashMap<String, usize>,
}

// rust only API
impl FragList {
    pub fn values(&self) -> impl Iterator<Item = &Fragment> {
        self.frag_list.values()
    }

    pub fn iter(&self) -> impl Iterator<Item = (&usize, &Fragment)> {
        self.frag_list.iter()
    }

    pub fn get(&self, frag_id: usize) -> Option<&Fragment> {
        self.frag_list.get(&frag_id)
    }
}

// Python+rust API
#[gen_stub_pymethods]
#[pymethods]
impl FragList {
    /// Create a new FragList, for n_atoms total atoms.
    ///
    /// Will create an ordered hashmap for fragments, and a dense Vec for atom_id -> fragment relationships.
    #[new]
    pub fn new(n_atoms: usize) -> Self {
        let defrag_list = (0..n_atoms).map(|_| Vec::new()).collect();
        Self {
            defrag_list,
            frag_list: OrderMap::new(),
            next_frag_id: 0,
            frag_counts: HashMap::new(), // used for rate calculation
        }
    }

    /// Add a fragment to frag list.
    ///
    /// Atom indices of -1 for missing optional and forbidden graph nodes.
    /// Will also update the defrag list.
    pub fn add_fragment(&mut self, name: String, atoms: Vec<isize>) -> usize {
        // frag counts
        if let Some(count) = self.frag_counts.get_mut(&name) {
            *count += 1;
        } else {
            self.frag_counts.insert(name.clone(), 1);
        };

        let id = self.next_frag_id;
        self.next_frag_id += 1;

        // defrag list
        for atom in atoms.iter() {
            if *atom == -1 {
                continue;
            }
            self.defrag_list.index_mut(*atom as usize).push(id);
        }

        // frag list
        let frag = Fragment {
            name,
            frag_id: id,
            atoms,
        };
        self.frag_list.insert(id, frag);
        id
    }

    /// Get a fragment corresponding to frag_id of index.
    ///
    /// Returns None if not found.
    pub fn get_fragment(&mut self, index: usize) -> Option<Fragment> {
        self.frag_list.get(&index).cloned()
    }

    /// Remove a fragment corresponding to frag_id of index.
    ///
    /// Returns whether a fragment of that frag_id was found.
    pub fn delete_fragment(&mut self, index: usize) -> bool {
        // frag list
        match self.frag_list.remove(&index) {
            Some(frag) => {
                // defrag list
                for atom in frag.atoms {
                    if atom == -1 {
                        continue;
                    }
                    self.defrag_list[atom as usize].retain(|e| *e != frag.frag_id)
                }
                // frag counts
                *self
                    .frag_counts
                    .get_mut(&frag.name)
                    .expect("delete_fragment but frag_counts did not contain key.") -= 1;
                true
            }
            None => false,
        }
    }

    /// Remove all fragments containing atoms.
    ///
    /// Called before recalculation of graphs for those atoms.
    pub fn delete_fragments_for_atoms(&mut self, atoms: Vec<usize>) {
        for atom in atoms {
            for frag_id in self.defrag_list[atom].clone() {
                self.delete_fragment(frag_id);
            }
        }
    }

    /// Get the number of fragments.
    ///
    /// Note: do not iterate 0 to this number, as frag_ids are not continuous.
    pub fn num_fragments(&self) -> usize {
        self.frag_list.len()
    }

    /// Get the number of fragments with a specific type.
    ///
    /// Always caches frag_counts, so should be fast.
    pub fn num_fragments_of_type(&self, name: &str) -> usize {
        self.frag_counts.get(name).map_or(0, |c| *c)
    }

    /// Return all frag_ids for a specific atom.
    ///
    /// Note: Will allocate a new list, but these lists are usually short, so should not be expensive.
    pub fn frag_ids_for(&self, atom_index: usize) -> Vec<usize> {
        self.defrag_list[atom_index].clone()
    }

    /// Get a copy of all frag IDs for debug purposes.
    ///
    /// Will make a copy, so should not be called often.
    pub fn get_all_frag_ids(&self) -> Vec<usize> {
        self.frag_list.keys().copied().collect()
    }
}
