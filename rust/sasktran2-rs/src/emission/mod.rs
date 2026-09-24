//! Band spectroscopy and population models, independent of VER interpolation
//! and radiative transfer. Chemistry-to-VER conversion belongs upstream.

pub mod o2;

use crate::atmosphere::StorageInputs;
use crate::prelude::*;

/// Fixed spectral data for one independently normalized photon emission band.
/// Level energies, transition probabilities, and excitation parameters belong
/// to the population model; the source only needs wavelengths and emitter mass.
#[derive(Clone, Debug)]
pub struct BandLineData {
    name: String,
    wavelengths_nm: Array1<f64>,
    molecular_mass_amu: f64,
}

impl BandLineData {
    pub fn new(
        name: impl Into<String>,
        wavelengths_nm: Array1<f64>,
        molecular_mass_amu: f64,
    ) -> Result<Self> {
        anyhow::ensure!(
            !wavelengths_nm.is_empty(),
            "Emission band must contain lines"
        );
        anyhow::ensure!(
            wavelengths_nm.iter().all(|w| w.is_finite() && *w > 0.0),
            "Emission wavelengths must be positive and finite"
        );
        anyhow::ensure!(
            molecular_mass_amu.is_finite() && molecular_mass_amu > 0.0,
            "Emission molecular mass must be positive and finite"
        );
        Ok(Self {
            name: name.into(),
            wavelengths_nm,
            molecular_mass_amu,
        })
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn wavelengths_nm(&self) -> ArrayView1<'_, f64> {
        self.wavelengths_nm.view()
    }

    pub fn molecular_mass_amu(&self) -> f64 {
        self.molecular_mass_amu
    }
}

/// Relative photon intensities within a band at the current atmospheric state.
///
/// Each weight array has shape (atmospheric location, line), with finite,
/// non-negative rows summing to one. The model owns the spectroscopy needed to
/// determine those weights and any excitation parameters or population profiles.
/// It need not assume rotational LTE or identify an excitation temperature with
/// the atmospheric kinetic temperature.
///
/// The temperature derivative is the local partial with respect to
/// `inputs.temperature_k()`, at fixed total band VER and fixed independent model
/// parameters. Its rows sum to zero. It excludes Doppler broadening, which is
/// differentiated by the common source using the kinetic temperature and mass.
/// A temperature-independent population model must explicitly return zero weight
/// derivatives; an omitted implementation must not silently imply zero.
///
/// The value-only method is used when temperature derivatives are disabled, so
/// implementations must not allocate or calculate derivatives on that path.
/// Static dispatch keeps native models out of Python and avoids virtual calls in
/// the line-profile loop.
pub trait BandEmissionModel {
    fn line_data(&self) -> &BandLineData;

    fn line_weights(&self, inputs: &impl StorageInputs) -> Result<Array2<f64>>;

    fn line_weights_with_temperature_derivative(
        &self,
        inputs: &impl StorageInputs,
    ) -> Result<(Array2<f64>, Array2<f64>)>;
}
