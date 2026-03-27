use ordermap::OrderMap;
use std::collections::HashMap;
use std::iter::Iterator;
use std::ops::IndexMut;
use pyo3::prelude::*;
use crate::fragment::Fragment;

#[pyclass]
pub struct FragList {
    n_atoms: usize,
    defrag_list: Vec<Vec<usize>>,
    frag_list: OrderMap<usize, Fragment>,
    #[pyo3(get)]
    next_frag_id: usize,
    #[pyo3(get)]
    frag_counts: HashMap<String, usize>,
}

// rust only API
impl FragList {
    pub fn values(&self) -> impl Iterator<Item = &Fragment>
    {
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
#[pymethods]
impl FragList {
    #[new]
    pub fn new(n_atoms: usize) -> Self {
        let defrag_list = (0..n_atoms).map(|_| Vec::new()).collect();
        Self {
            n_atoms,
            defrag_list,
            frag_list: OrderMap::new(),
            next_frag_id: 0,
            frag_counts: HashMap::new(), // used for rate calculation
        }
    }

    pub fn add_fragment(
        &mut self,
        name: String,
        atoms: Vec<isize>,
    ) -> usize {
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

    pub fn get_fragment(&mut self, index: usize) -> Option<Fragment> {
        let frag = self.frag_list.get(&index);
        match frag {
            // returns a read only copy
            // usually there isn't that many atoms per fragment, so this isn't so bad
            Some(frag) => Some(frag.clone()),
            None => None,
        }
    }

    pub fn delete_fragment(&mut self, index: usize) -> bool {
        // frag list
        match self.frag_list.remove(&index) {
            Some(frag) => {
                // defrag list
                for atom in frag.atoms {
                    if atom == -1 {
                        continue
                    }
                    self.defrag_list[atom as usize].retain(
                        |e| *e != frag.frag_id
                    )
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

    pub fn delete_fragments_for_atoms(&mut self, atoms: Vec<usize>) {
        for atom in atoms {
            for frag_id in self.defrag_list[atom].clone() {
                self.delete_fragment(frag_id);
            }
        }
    }

    pub fn num_fragments(&self) -> usize {
        self.frag_list.len()
    }

    pub fn num_fragments_of_type(&self, name: &str) -> usize {
        self.frag_counts.get(name).map_or(0, |c| *c)
    }

    pub fn frag_ids_for(&self, atom_index: usize) -> Vec<usize> {
        // there usually aren't that many fragments per atom, so this isn't so bad
        self.defrag_list[atom_index].clone()
    }
}
