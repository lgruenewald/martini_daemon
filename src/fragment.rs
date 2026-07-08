use std::fmt::Display;

use pyo3::prelude::*;
use pyo3_stub_gen::derive::gen_stub_pyclass;

#[gen_stub_pyclass]
#[pyclass(str, skip_from_py_object)]
#[derive(Clone)]
pub struct Fragment {
    // must be same as graph name
    #[pyo3(get)]
    pub name: String,
    #[pyo3(get)]
    pub frag_id: usize,
    // missing optional -> -1
    #[pyo3(get)]
    pub atoms: Vec<isize>,
}

#[pymethods]
impl Fragment {
    #[new]
    pub fn new(name: String, frag_id: usize, atoms: Vec<isize>) -> Self {
        Self {
            name, frag_id, atoms
        }
    }
}

impl Display for Fragment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}({}; {:?})", self.name, self.frag_id, self.atoms)
    }
}
