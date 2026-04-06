use pyo3::prelude::*;
use std::env;

mod periodic_box;
mod fragment;
mod detection_template;
mod detection;
mod detection_one;
mod frag_list;
mod detection_template_list;
mod bond_graph;

#[pyfunction]
fn build_version() -> String {
    env!("BUILD_VERSION").to_string()
}

#[pymodule]
#[pyo3(name = "__rust")]
fn rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<periodic_box::PeriodicBox>()?;
    m.add_class::<fragment::Fragment>()?;
    m.add_class::<frag_list::FragList>()?;
    m.add_class::<detection_template::DetectionTemplate>()?;
    m.add_class::<detection_template_list::DetectionTemplateList>()?;
    m.add_class::<bond_graph::BondGraph>()?;
    m.add_function(wrap_pyfunction!(detection::detection, m)?)?;
    m.add_function(wrap_pyfunction!(build_version, m)?)?;
    Ok(())
}
