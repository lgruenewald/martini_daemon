use pyo3::prelude::*;
use pyo3_stub_gen::{define_stub_info_gatherer, derive::gen_stub_pyfunction};
use std::env;

mod bond_graph;
mod detection;
mod detection_one;
mod detection_template;
mod detection_template_list;
mod frag_list;
mod fragment;
mod parser;
mod periodic_box;
mod toptraj;

#[gen_stub_pyfunction]
#[pyfunction]
/// Return the current version and git commit as a string.
///
/// The output gets embedded in logs.
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
    m.add_class::<toptraj::TopTrajFrame>()?;
    m.add_class::<toptraj::TopTrajReader>()?;
    m.add_class::<toptraj::TopTrajWriter>()?;
    m.add_function(wrap_pyfunction!(detection::detection, m)?)?;
    m.add_function(wrap_pyfunction!(build_version, m)?)?;
    m.add_function(wrap_pyfunction!(parser::tokenize, m)?)?;
    Ok(())
}

define_stub_info_gatherer!(stub_info);
