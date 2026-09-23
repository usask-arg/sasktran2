use numpy::*;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;
use std::path::PathBuf;
use std::str::FromStr;

use sasktran2_rs::constituent::traits::Constituent;
use sasktran2_rs::constituent::types::band_volume_emission_rate::{
    BandVolumeEmissionRate, oxygen_emission_band, validate_photon_ver,
};
use sasktran2_rs::optical::line::hitran_loader::{hitran_molecule_file, read_hitran_line_file};
use sasktran2_rs::photchem::emission::AEmissionLineWeightModel;

use crate::constituent::atmo_storage::AtmosphereStorage;

#[pyclass]
pub struct PyO2BandEmissionRate {
    pub inner: BandVolumeEmissionRate,
}

#[pymethods]
impl PyO2BandEmissionRate {
    #[new]
    #[pyo3(
        signature = (altitudes_m, photon_ver, hitran_directory, band = "0-0", line_weight_model = "einstein_a_branching", out_of_bounds_mode = "zero"),
    )]
    fn new<'py>(
        altitudes_m: PyReadonlyArray1<'py, f64>,
        photon_ver: PyReadonlyArray1<'py, f64>,
        hitran_directory: &str,
        band: &str,
        line_weight_model: &str,
        out_of_bounds_mode: Option<&str>,
    ) -> PyResult<Self> {
        let line_weight_model = AEmissionLineWeightModel::from_str(line_weight_model)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let path = hitran_molecule_file("O2", &PathBuf::from(hitran_directory))
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let db = read_hitran_line_file(path).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let band =
            oxygen_emission_band(&db, band).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let mut inner = BandVolumeEmissionRate::new(
            altitudes_m.as_array().to_owned(),
            photon_ver.as_array().to_owned(),
            band,
            line_weight_model,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

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
            };
        }

        Ok(Self { inner })
    }

    #[getter]
    fn get_photon_ver<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        let array = &this.borrow().inner.photon_ver;

        unsafe { PyArray1::borrow_from_array(array, this.into_any()) }
    }

    #[setter]
    fn set_photon_ver(&mut self, photon_ver: PyReadonlyArray1<f64>) -> PyResult<()> {
        // The input may itself be a view of our current profile.
        let photon_ver = photon_ver.as_array().to_owned();
        validate_photon_ver(photon_ver.view(), self.inner.altitudes.len())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        // Preserve any previously returned mutable numpy view of this profile.
        self.inner.photon_ver.assign(&photon_ver);

        Ok(())
    }

    #[getter]
    fn get_altitudes_m<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.altitudes.clone().into_pyarray(py)
    }

    #[getter]
    fn get_wavelengths_nm<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.band.wavelengths_nm().into_pyarray(py)
    }

    #[getter]
    fn get_band(&self) -> &str {
        &self.inner.band.name
    }

    fn line_weights<'py>(
        &self,
        py: Python<'py>,
        temperature_k: PyReadonlyArray1<f64>,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        self.inner
            .line_weights(temperature_k.as_array())
            .map(|weights| weights.into_pyarray(py))
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    pub fn add_to_atmosphere<'py>(&mut self, atmo: Bound<'py, PyAny>) -> PyResult<()> {
        let mut rust_atmo = AtmosphereStorage::new(&atmo)?;

        self.inner
            .add_to_atmosphere(&mut rust_atmo)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

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
