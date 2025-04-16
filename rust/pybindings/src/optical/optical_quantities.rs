use ndarray::*;
use numpy::*;
use pyo3::prelude::*;
use sk_core::optical::OpticalQuantities;

#[pyclass]
pub struct PyOpticalQuantities {
    oq: OpticalQuantities,
}

impl PyOpticalQuantities {
    pub fn new(oq: OpticalQuantities) -> Self {
        Self { oq }
    }
}

#[pymethods]
impl PyOpticalQuantities {
    #[new]
    fn py_new(num_geometry: usize, num_wavelengths: usize) -> Self {
        let oq = OpticalQuantities::new(num_geometry, num_wavelengths);
        Self { oq }
    }

    #[getter]
    fn get_cross_section<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray2<f64>> {
        let array = &this.borrow().oq.cross_section;

        unsafe { PyArray2::borrow_from_array(array, this.into_any()) }
    }

    #[getter]
    // This is just for interfacing
    fn get_extinction<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray2<f64>> {
        let array = &this.borrow().oq.cross_section;

        unsafe { PyArray2::borrow_from_array(array, this.into_any()) }
    }

    #[getter]
    fn get_ssa<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray2<f64>> {
        let array = &this.borrow().oq.ssa;

        unsafe { PyArray2::borrow_from_array(array, this.into_any()) }
    }
}
