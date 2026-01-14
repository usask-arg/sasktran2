use pyo3::prelude::*;
use sasktran2_rs::bindings::common;

#[pyfunction]
pub fn openmp_support_enabled() -> PyResult<bool> {
    Ok(common::openmp_support_enabled())
}

#[pyfunction]
pub fn lto_test() -> PyResult<bool> {
    common::lto_test();
    Ok(true)
}
