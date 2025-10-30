use numpy::{PyArray2, PyArray3, PyReadonlyArray2, PyReadonlyArray3};
use pyo3::prelude::*;

use sasktran2_rs::constituent::traits::Constituent;
use sasktran2_rs::constituent::types::manual::Manual as RustManual;

use crate::constituent::atmo_storage::AtmosphereStorage;

#[pyclass]
pub struct PyManual {
    inner: RustManual,
}

#[pymethods]
impl PyManual {
    #[new]
    fn new<'py>(
        extinction: PyReadonlyArray2<'py, f64>,
        ssa: PyReadonlyArray2<'py, f64>,
        legendre_moments: Option<PyReadonlyArray3<'py, f64>>,
    ) -> Self {
        if let Some(legendre_moments) = legendre_moments {
            Self {
                inner: RustManual::new(
                    extinction.as_array().to_owned(),
                    ssa.as_array().to_owned(),
                    Some(legendre_moments.as_array().to_owned()),
                ),
            }
        } else {
            Self {
                inner: RustManual::new(
                    extinction.as_array().to_owned(),
                    ssa.as_array().to_owned(),
                    None,
                ),
            }
        }
    }

    pub fn add_to_atmosphere<'py>(&self, atmo: &Bound<'py, PyAny>) -> PyResult<()> {
        let mut rust_atmo = AtmosphereStorage::new(atmo)?;

        let _ = self.inner.add_to_atmosphere(&mut rust_atmo);

        Ok(())
    }

    pub fn register_derivative(&self, atmo: &'_ Bound<'_, PyAny>, name: &str) -> PyResult<()> {
        let mut rust_atmo = AtmosphereStorage::new(atmo)?;

        let _ = self.inner.register_derivatives(&mut rust_atmo, name);

        Ok(())
    }

    #[getter]
    fn get_extinction<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray2<f64>> {
        let array = &this.borrow().inner.extinction;

        unsafe { PyArray2::borrow_from_array(array, this.into_any()) }
    }

    #[getter]
    fn get_ssa<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray2<f64>> {
        let array = &this.borrow().inner.ssa;

        unsafe { PyArray2::borrow_from_array(array, this.into_any()) }
    }

    #[getter]
    fn get_leg_coeff<'py>(this: Bound<'py, Self>) -> Option<Bound<'py, PyArray3<f64>>> {
        if let Some(legendre_moments) = &this.borrow().inner.legendre_moments {
            let array = legendre_moments;
            unsafe { Some(PyArray3::borrow_from_array(array, this.into_any())) }
        } else {
            None
        }
    }
}
