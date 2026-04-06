use std::fmt::Display;

use pyo3::prelude::*;

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

impl Display for Fragment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}({}; {:?})", self.name, self.frag_id, self.atoms)
    }
}
