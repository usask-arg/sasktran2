//! Rayleigh scattering parametirizations and interfaces
//! to SASKTRAN2.  This is primarly the Bates parameterization for the
//! cross section and King factor, and then a Constituent interface to SASKTRAN2

use numpy::*;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;

use sk_core::constituent::traits::Constituent;
use sk_core::constituent::types::vmr_alt_absorber::VMRAltitudeAbsorber as VMRAltitudeAbsorberCore;

use crate::constituent::atmo_storage::AtmosphereStorage;
use crate::optical::optical_property::PyOpticalProperty;

#[pyclass]
pub struct PyVMRAltitudeAbsorber {
    pub inner: VMRAltitudeAbsorberCore<PyOpticalProperty>,
    optical_property: Py<PyAny>,
}

#[pymethods]
impl PyVMRAltitudeAbsorber {
    #[new]
    #[pyo3(
        signature = (optical_property, altitudes_m, vmr, out_of_bounds_mode = "zero"),
    )]
    fn new<'py>(
        optical_property: Bound<'py, PyAny>,
        altitudes_m: PyReadonlyArray1<'py, f64>,
        vmr: PyReadonlyArray1<'py, f64>,
        out_of_bounds_mode: Option<&str>,
    ) -> PyResult<Self> {
        let mut inner = VMRAltitudeAbsorberCore::new(
            altitudes_m.as_array().to_owned(),
            vmr.as_array().to_owned(),
        );

        if out_of_bounds_mode.is_some() {
            inner = match out_of_bounds_mode.unwrap() {
                "zero" => inner.with_interp_mode(sk_core::interpolation::OutOfBoundsMode::Zero),
                "extend" => inner.with_interp_mode(sk_core::interpolation::OutOfBoundsMode::Extend),
                mode => {
                    return Err(PyValueError::new_err(format!(
                        "Invalid out_of_bounds_mode: {}",
                        mode
                    )));
                }
            }
        }

        Ok(Self {
            inner: inner,
            optical_property: optical_property.into(),
        })
    }

    #[getter]
    fn get_vmr<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        let array = &this.borrow().inner.vmr;

        unsafe { PyArray1::borrow_from_array(array, this.into_any()) }
    }

    #[getter]
    fn get_altitudes_m<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        let array = &this.borrow().inner.altitudes;

        unsafe { PyArray1::borrow_from_array(array, this.into_any()) }
    }

    ///
    /// Test docstring add to atmosphere
    pub fn add_to_atmosphere<'py>(&mut self, atmo: Bound<'py, PyAny>) -> PyResult<()> {
        let mut rust_atmo = AtmosphereStorage::new(&atmo);

        let py_optical =
            PyOpticalProperty::new(self.optical_property.clone_ref(atmo.py()), atmo.into());

        let _ = self.inner.with_optical_property(py_optical);
        let _ = self.inner.add_to_atmosphere(&mut rust_atmo);

        Ok(())
    }

    pub fn register_derivative(&mut self, atmo: Bound<'_, PyAny>, name: &str) -> PyResult<()> {
        let mut rust_atmo = AtmosphereStorage::new(&atmo);

        let py_optical =
            PyOpticalProperty::new(self.optical_property.clone_ref(atmo.py()), atmo.into());

        let _ = self.inner.with_optical_property(py_optical);

        self.inner
            .register_derivatives(&mut rust_atmo, name)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(())
    }
}
