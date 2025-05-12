use pyo3::prelude::*;
use sasktran2_rs::bindings::common;

#[pyfunction]
pub fn openmp_support_enabled() -> PyResult<bool> {
    Ok(common::openmp_support_enabled())
}
