//! Emission

use pyo3::prelude::*;
use pyo3::types::PyAny;

use sasktran2_rs::constituent::traits::Constituent;
use sasktran2_rs::constituent::types::emission::ThermalEmission as RustThermalEmissionCore;

use crate::constituent::atmo_storage::AtmosphereStorage;

#[pyclass]
/// An implementation of thermal emissions calculated from the Planck function. The emission is
/// calculated with units of [W / (m^2 nm sr)].
///
/// This Constituent requires that the atmosphere object have `temperature_k` and
/// `wavelength_nm` defined inside the :py:class:`sasktran2.Atmosphere` object.
pub struct PyThermalEmission {
    inner: RustThermalEmissionCore,
}

#[pymethods]
impl PyThermalEmission {
    #[new]
    #[pyo3(
        signature = (),
    )]
    fn new() -> Self {
        let inner = RustThermalEmissionCore::new();

        PyThermalEmission { inner }
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
}
