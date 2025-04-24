use pyo3::prelude::*;
use sasktran2_rs::bindings::deriv_mapping;
use numpy::*;

#[pyclass(unsendable)]
pub struct PyDerivativeMappingView {
    pub derivative_mapping: deriv_mapping::DerivativeMapping,
}

#[pymethods]
impl PyDerivativeMappingView {
    #[getter]
    fn get_d_ssa<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let binding = this.borrow();
        let array = &binding.derivative_mapping.d_ssa();

        unsafe { Ok(PyArray2::borrow_from_array(array, this.into_any())) }
    }
}