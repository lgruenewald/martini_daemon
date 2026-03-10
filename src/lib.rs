use pyo3::prelude::*;

mod periodic_box;

#[pymodule]
#[pyo3(name = "__rust")]
fn rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<periodic_box::PeriodicBox>()?;
    Ok(())
}

