use numpy::PyArray3;
use pyo3::prelude::*;
use sasktran2_rs::bindings::output;

#[pyclass(unsendable)]
pub struct PyOutput {
    pub output: output::Output,
}

impl PyOutput {
    pub fn new(num_wavel: usize, num_los: usize, num_stokes: usize) -> Self {
        let output = output::Output::new(num_wavel, num_los, num_stokes);
        Self { output }
    }
}

#[pymethods]
impl PyOutput {
    #[getter]
    fn get_radiance<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray3<f64>>> {
        let array = &this.borrow().output.radiance;

        unsafe { Ok(PyArray3::borrow_from_array(array, this.into_any())) }
    }
}