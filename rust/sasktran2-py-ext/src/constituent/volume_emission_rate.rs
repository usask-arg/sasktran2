//! Rayleigh scattering parametirizations and interfaces
//! to SASKTRAN2.  This is primarly the Bates parameterization for the
//! cross section and King factor, and then a Constituent interface to SASKTRAN2

use numpy::*;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;

use sasktran2_rs::constituent::traits::Constituent;
use sasktran2_rs::constituent::types::volume_emission_rate::MonochromaticVolumeEmissionRate;

use crate::constituent::atmo_storage::AtmosphereStorage;

#[pyclass]
pub struct PyMonochromaticVolumeEmissionRate {
    pub inner: MonochromaticVolumeEmissionRate,
}

#[pymethods]
impl PyMonochromaticVolumeEmissionRate {
    #[new]
    #[pyo3(
        signature = (altitudes_m, ver, wavelength_nm, out_of_bounds_mode = "zero"),
    )]
    fn new<'py>(
        altitudes_m: PyReadonlyArray1<'py, f64>,
        ver: PyReadonlyArray1<'py, f64>,
        wavelength_nm: f64,
        out_of_bounds_mode: Option<&str>,
    ) -> PyResult<Self> {
        let mut inner = MonochromaticVolumeEmissionRate::new(
            altitudes_m.as_array().to_owned(),
            ver.as_array().to_owned(),
            wavelength_nm,
        );

        if let Some(out_of_bounds_mode) = out_of_bounds_mode {
            inner = match out_of_bounds_mode {
                "zero" => {
                    inner.with_interp_mode(sasktran2_rs::interpolation::OutOfBoundsMode::Zero)
                }
                "extend" => {
                    inner.with_interp_mode(sasktran2_rs::interpolation::OutOfBoundsMode::Extend)
                }
                mode => {
                    return Err(PyValueError::new_err(format!(
                        "Invalid out_of_bounds_mode: {mode}"
                    )));
                }
            }
        }

        Ok(Self { inner })
    }

    #[getter]
    fn get_ver<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        let array = &this.borrow().inner.ver;

        unsafe { PyArray1::borrow_from_array(array, this.into_any()) }
    }

    #[setter]
    fn set_ver(&mut self, ver: PyReadonlyArray1<f64>) -> PyResult<()> {
        // Check if the length of the new vmr matches the length of the altitudes
        if ver.len() != self.inner.altitudes.len() {
            return Err(PyValueError::new_err(format!(
                "Length of vmr ({}) does not match length of altitudes ({})",
                ver.len(),
                self.inner.altitudes.len()
            )));
        }

        self.inner.ver = ver.as_array().to_owned();

        Ok(())
    }

    #[getter]
    fn get_altitudes_m<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        let array = &this.borrow().inner.altitudes;

        unsafe { PyArray1::borrow_from_array(array, this.into_any()) }
    }

    ///
    /// Test docstring add to atmosphere
    pub fn add_to_atmosphere<'py>(&mut self, atmo: Bound<'py, PyAny>) -> PyResult<()> {
        let mut rust_atmo = AtmosphereStorage::new(&atmo)?;

        let _ = self.inner.add_to_atmosphere(&mut rust_atmo);

        Ok(())
    }

    pub fn register_derivative(&mut self, atmo: Bound<'_, PyAny>, name: &str) -> PyResult<()> {
        let mut rust_atmo = AtmosphereStorage::new(&atmo)?;

        self.inner
            .register_derivatives(&mut rust_atmo, name)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(())
    }
}
