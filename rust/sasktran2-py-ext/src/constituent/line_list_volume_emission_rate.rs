use ndarray::{Array2, Ix1, Ix2};
use numpy::*;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;

use sasktran2_rs::constituent::traits::Constituent;
use sasktran2_rs::constituent::types::line_list_volume_emission_rate::LineListVolumeEmissionRate;

use crate::constituent::atmo_storage::AtmosphereStorage;

#[pyclass]
pub struct PyLineListVolumeEmissionRate {
    pub inner: LineListVolumeEmissionRate,
}

#[pymethods]
impl PyLineListVolumeEmissionRate {
    #[new]
    #[pyo3(
        signature = (altitudes_m, photon_ver, wavelengths_nm, weights, out_of_bounds_mode = "zero"),
    )]
    fn new<'py>(
        altitudes_m: PyReadonlyArray1<'py, f64>,
        photon_ver: PyReadonlyArray1<'py, f64>,
        wavelengths_nm: PyReadonlyArray1<'py, f64>,
        weights: PyReadonlyArrayDyn<'py, f64>,
        out_of_bounds_mode: Option<&str>,
    ) -> PyResult<Self> {
        let weights = weights_from_py(weights, altitudes_m.len())?;

        let mut inner = LineListVolumeEmissionRate::new(
            altitudes_m.as_array().to_owned(),
            photon_ver.as_array().to_owned(),
            wavelengths_nm.as_array().to_owned(),
            weights,
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
        if photon_ver.len() != self.inner.altitudes.len() {
            return Err(PyValueError::new_err(format!(
                "Length of photon_ver ({}) does not match length of altitudes ({})",
                photon_ver.len(),
                self.inner.altitudes.len()
            )));
        }

        self.inner.photon_ver = photon_ver.as_array().to_owned();

        Ok(())
    }

    #[getter]
    fn get_altitudes_m<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        let array = &this.borrow().inner.altitudes;

        unsafe { PyArray1::borrow_from_array(array, this.into_any()) }
    }

    #[getter]
    fn get_wavelengths_nm<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        let array = &this.borrow().inner.wavelengths_nm;

        unsafe { PyArray1::borrow_from_array(array, this.into_any()) }
    }

    #[getter]
    fn get_weights<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray2<f64>> {
        let array = &this.borrow().inner.weights;

        unsafe { PyArray2::borrow_from_array(array, this.into_any()) }
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

fn weights_from_py(
    weights: PyReadonlyArrayDyn<'_, f64>,
    num_altitude: usize,
) -> PyResult<Array2<f64>> {
    match weights.ndim() {
        1 => {
            let weights_1d = weights
                .as_array()
                .into_dimensionality::<Ix1>()
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
            Ok(Array2::from_shape_fn(
                (num_altitude, weights_1d.len()),
                |(_, line_idx)| weights_1d[line_idx],
            ))
        }
        2 => weights
            .as_array()
            .into_dimensionality::<Ix2>()
            .map(|weights| weights.to_owned())
            .map_err(|e| PyValueError::new_err(e.to_string())),
        ndim => Err(PyValueError::new_err(format!(
            "weights must be one- or two-dimensional, got {ndim} dimensions"
        ))),
    }
}
