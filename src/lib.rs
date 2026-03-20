use pyo3::prelude::*;
use std::env;

mod periodic_box;

#[pyfunction]
fn build_version() -> String {
    env!("BUILD_VERSION").to_string()
}

#[pymodule]
#[pyo3(name = "__rust")]
fn rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<periodic_box::PeriodicBox>()?;
    m.add_function(wrap_pyfunction!(build_version, m)?)?;
    Ok(())
}

