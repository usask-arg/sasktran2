//! Rayleigh scattering parametirizations and interfaces
//! to SASKTRAN2.  This is primarly the Bates parameterization for the
//! cross section and King factor, and then a Constituent interface to SASKTRAN2

use numpy::PyReadonlyArray1;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;

use sk_core::constituent::rayleigh::Rayleigh as RustRayleighCore;
use sk_core::constituent::{Constituent, StorageInputs};

use crate::constituent::atmo_storage::AtmosphereStorage;

#[pyclass]
/// An implementation of Rayleigh scattering.  Cross sections (and depolarization factors) can be
/// calculated multiple ways, with the default method being that of 'bates'.
///
/// Rayleigh scattering number density is estimated through the ideal gas law.
///
/// This Constituent requires that the atmosphere object have `temperature_k`, `pressure_pa`, and
/// `wavelength_nm` are all defined inside the :py:class:`sasktran2.Atmosphere` object.
///
/// Parameters
/// ----------
/// method : str, default='bates'
///     Method to use to calculate the cross section.  Supported methods are
///     ['bates', 'manual'], by default 'bates'
/// n2_percentage : float, optional
///    Percentage of N2 in the atmosphere, by default 78.084
/// o2_percentage : float, optional
///     Percentage of O2 in the atmosphere, by default 20.946
/// ar_percentage : float, optional
///    Percentage of Ar in the atmosphere, by default 0.934
/// co2_percentage : float, optional
///    Percentage of CO2 in the atmosphere, by default 0.036
/// wavelengths_nm : numpy.ndarray, optional
///    Wavelengths in nm to use for the cross section, by default None, only used when method is "manual"
/// xs : numpy.ndarray, optional
///    Cross section in m2/molecule to use for the cross section, by default None, only used when method is "manual"
/// king : numpy.ndarray, optional
///    King factor to use for the cross section, by default None, only used when method is "manual"

pub struct Rayleigh {
    inner: RustRayleighCore,
}

#[pymethods]
impl Rayleigh {
    #[new]
    #[pyo3(
        signature = (method="bates", n2_percentage=None, o2_percentage=None, ar_percentage=None, co2_percentage=None, wavelengths_nm=None, xs=None, king=None),
    )]
    /// Test dosstring new
    fn new<'py>(
        method: &str,
        n2_percentage: Option<f64>,
        o2_percentage: Option<f64>,
        ar_percentage: Option<f64>,
        co2_percentage: Option<f64>,
        wavelengths_nm: Option<PyReadonlyArray1<'py, f64>>,
        xs: Option<PyReadonlyArray1<'py, f64>>,
        king: Option<PyReadonlyArray1<'py, f64>>,
    ) -> Self {
        let mut inner = RustRayleighCore::new();

        if method.to_lowercase() == "manual" {
            let xs = xs
                .ok_or_else(|| {
                    PyErr::new::<pyo3::exceptions::PyValueError, _>(
                        "xs must be specified when using the manual method",
                    )
                })
                .unwrap()
                .as_array()
                .to_owned();

            let king = king
                .ok_or_else(|| {
                    PyErr::new::<pyo3::exceptions::PyValueError, _>(
                        "king must be specified when using the manual method",
                    )
                })
                .unwrap()
                .as_array()
                .to_owned();

            let wavelengths_nm = wavelengths_nm
                .ok_or_else(|| {
                    PyErr::new::<pyo3::exceptions::PyValueError, _>(
                        "wavelengths_nm must be specified when using the manual method",
                    )
                })
                .unwrap()
                .as_array()
                .to_owned();

            inner = inner.with_manual_xs(xs, king, wavelengths_nm);
        }

        if n2_percentage.is_some() {
            inner = inner.with_n2_percentage(n2_percentage.unwrap());
        }

        if o2_percentage.is_some() {
            inner = inner.with_o2_percentage(o2_percentage.unwrap());
        }

        if ar_percentage.is_some() {
            inner = inner.with_ar_percentage(ar_percentage.unwrap());
        }

        if co2_percentage.is_some() {
            inner = inner.with_co2_percentage(co2_percentage.unwrap());
        }

        Rayleigh { inner: inner }
    }

    ///
    /// Test docstring add to atmosphere
    pub fn add_to_atmosphere<'py>(&self, atmo: &Bound<'py, PyAny>) -> PyResult<()> {
        let rust_atmo = AtmosphereStorage::new(atmo);

        let inputs = rust_atmo.inputs;
        let mut outputs = rust_atmo.outputs;

        inputs.wavelengths_nm().ok_or_else(|| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "wavelengths_nm must be specified in the atmosphere to use this function",
            )
        })?;

        if rust_atmo.num_stokes != 1 && rust_atmo.num_stokes != 3 {
            return Err(PyValueError::new_err(format!(
                "Invalid number of Stokes parameters: {} (expected 1 or 3)",
                rust_atmo.num_stokes
            )));
        }
        self.inner.add_to_atmosphere(&inputs, &mut outputs);

        Ok(())
    }

    pub fn register_derivative(&self, atmo: &'_ Bound<'_, PyAny>, name: &str) -> PyResult<()> {
        let rust_atmo = AtmosphereStorage::new(atmo);

        let inputs = rust_atmo.inputs;
        let outputs = rust_atmo.outputs;
        let deriv_generator = rust_atmo.deriv_generator;

        self.inner
            .register_derivatives(&inputs, &outputs, name, &deriv_generator);

        Ok(())
    }
}
