//! An independently specified band photon VER with a replaceable population
//! model. Spectroscopy and chemistry-to-VER conversion belong upstream.
use crate::atmosphere::types::LineshapeCoordinate;
use crate::atmosphere::*;
use crate::bindings::config::SpectralGridMode;
use crate::constituent::traits::Constituent;
use crate::emission::BandEmissionModel;
use crate::interpolation::{OutOfBoundsMode, linear::linear_interpolating_matrix};
use crate::math::errorfunctions::optimized::SQRT_PI;
use crate::optical::line::shape::{LineShape, LineShapeDirection, assign_with_derivative};
use crate::optical::line::{AdjustedLineParameters, OpticalLine};
use crate::optical::types::line_absorber::assign_normalized_doppler_line_shape;
use crate::prelude::*;

const FOUR_PI: f64 = 4.0 * std::f64::consts::PI;

/// One band per constituent gives a separate VER retrieval parameter, even
/// when several vibrational bands overlap spectrally.
#[derive(Clone)]
pub struct BandVolumeEmissionRate<M: BandEmissionModel> {
    pub altitudes: Array1<f64>,
    pub photon_ver: Array1<f64>,
    pub model: M,
    interp_mode: OutOfBoundsMode,
}

impl<M: BandEmissionModel> BandVolumeEmissionRate<M> {
    pub fn new(altitudes: Array1<f64>, photon_ver: Array1<f64>, model: M) -> Result<Self> {
        anyhow::ensure!(
            !altitudes.is_empty(),
            "Emission altitude grid cannot be empty"
        );
        anyhow::ensure!(
            altitudes.iter().all(|a| a.is_finite())
                && altitudes
                    .iter()
                    .zip(altitudes.iter().skip(1))
                    .all(|(a, b)| b > a),
            "Emission altitudes must be finite and strictly increasing"
        );
        validate_photon_ver(photon_ver.view(), altitudes.len())?;
        Ok(Self {
            altitudes,
            photon_ver,
            model,
            interp_mode: OutOfBoundsMode::Zero,
        })
    }

    pub fn with_interp_mode(mut self, mode: OutOfBoundsMode) -> Self {
        self.interp_mode = mode;
        self
    }

    pub fn band_name(&self) -> &str {
        self.model.line_data().name()
    }

    pub fn wavelengths_nm(&self) -> ArrayView1<'_, f64> {
        self.model.line_data().wavelengths_nm()
    }

    fn interpolation(&self, inputs: &impl StorageInputs) -> Array2<f64> {
        linear_interpolating_matrix(&self.altitudes, &inputs.altitude_m(), self.interp_mode)
    }

    fn spectrum(
        &self,
        inputs: &impl StorageInputs,
        derivative: bool,
        ver: Option<ArrayView1<'_, f64>>,
    ) -> Result<(Array2<f64>, Option<Array2<f64>>)> {
        let (mut weights, mut d_weights) = if derivative {
            let (weights, derivatives) = self
                .model
                .line_weights_with_temperature_derivative(inputs)?;
            (weights, Some(derivatives))
        } else {
            (self.model.line_weights(inputs)?, None)
        };
        let lines = self.model.line_data();
        let expected = (inputs.altitude_m().len(), lines.wavelengths_nm().len());
        anyhow::ensure!(
            weights.dim() == expected,
            "Band line weights must have shape (atmospheric locations, lines)"
        );
        anyhow::ensure!(
            d_weights.as_ref().is_none_or(|d| d.dim() == expected),
            "Band line-weight temperature derivatives must match the weight shape"
        );
        if let Some(ver) = ver {
            Zip::from(weights.rows_mut())
                .and(ver)
                .for_each(|mut row, &v| row *= v);
            if let Some(d_weights) = d_weights.as_mut() {
                Zip::from(d_weights.rows_mut())
                    .and(ver)
                    .for_each(|mut row, &v| row *= v);
            }
        }
        doppler_spectrum(
            inputs,
            lines.wavelengths_nm(),
            lines.molecular_mass_amu(),
            weights.view(),
            d_weights.as_ref().map(|d| d.view()),
        )
    }

    /// Population wrappers request T at fixed populations/Einstein coefficients;
    /// direct VER constituents additionally expose their own VER mapping.
    pub fn register_emission_derivatives(
        &self,
        storage: &mut impl AtmosphereStorageAccess,
        name: &str,
        register_ver: bool,
    ) -> Result<()> {
        let (inputs, _, generator) = storage.split_inputs_outputs_deriv();
        let temperature_derivative = inputs.calculate_temperature_derivative();
        if !register_ver && !temperature_derivative {
            return Ok(());
        }
        validate_photon_ver(self.photon_ver.view(), self.altitudes.len())?;
        let interpolation = self.interpolation(inputs);
        let ver = interpolation.dot(&self.photon_ver);
        let (spectrum, d_spectrum) = self.spectrum(
            inputs,
            temperature_derivative,
            if register_ver { None } else { Some(ver.view()) },
        )?;
        if register_ver {
            let full_name = format!("wf_{name}_photon_ver");
            let mut mapping = generator.get_derivative_mapping(&full_name);
            mapping.set_assign_name(&full_name);
            mapping.set_interp_dim(&format!("{name}_altitude"));
            mapping.set_interpolator(&interpolation);
            mapping.mut_view().d_emission.assign(&(&spectrum / FOUR_PI));
        }
        if let Some(mut d_spectrum) = d_spectrum {
            Zip::from(d_spectrum.rows_mut())
                .and(&ver)
                .for_each(|mut row, &v| {
                    row *= if register_ver {
                        v / FOUR_PI
                    } else {
                        1.0 / FOUR_PI
                    }
                });
            let full_name = format!("wf_{name}_temperature_k");
            let mut mapping = generator.get_derivative_mapping(&full_name);
            mapping.set_assign_name("wf_temperature_k");
            mapping.set_interp_dim("altitude");
            mapping.mut_view().d_emission.assign(&d_spectrum);
        }
        Ok(())
    }
}

impl<M: BandEmissionModel> Constituent for BandVolumeEmissionRate<M> {
    fn add_to_atmosphere(&self, storage: &mut impl AtmosphereStorageAccess) -> Result<()> {
        let (inputs, outputs) = storage.split_inputs_outputs();
        validate_photon_ver(self.photon_ver.view(), self.altitudes.len())?;
        let ver = self.interpolation(inputs).dot(&self.photon_ver);
        let (mut spectrum, _) = self.spectrum(inputs, false, Some(ver.view()))?;
        spectrum /= FOUR_PI;
        outputs.mut_view().emission_source += &spectrum;
        Ok(())
    }

    fn register_derivatives(
        &self,
        storage: &mut impl AtmosphereStorageAccess,
        name: &str,
    ) -> Result<()> {
        self.register_emission_derivatives(storage, name, true)
    }
}

pub fn validate_photon_ver(ver: ArrayView1<'_, f64>, size: usize) -> Result<()> {
    anyhow::ensure!(
        ver.len() == size,
        "Photon VER length must match emission altitudes"
    );
    anyhow::ensure!(
        ver.iter().all(|v| v.is_finite() && *v >= 0.0),
        "Photon VER must be finite and non-negative"
    );
    Ok(())
}

/// Values use the original Gaussian kernel; derivative arrays and directional
/// kernels are only used when d_areas is supplied. A zero d_areas includes just
/// Doppler broadening, whereas band weights also contribute an amplitude term.
pub(crate) fn doppler_spectrum(
    inputs: &impl StorageInputs,
    wavelengths: ArrayView1<'_, f64>,
    molecular_mass_amu: f64,
    areas: ArrayView2<'_, f64>,
    d_areas: Option<ArrayView2<'_, f64>>,
) -> Result<(Array2<f64>, Option<Array2<f64>>)> {
    let temperature = inputs
        .temperature_k()
        .ok_or_else(|| anyhow!("Temperature must be set for line emission"))?;
    anyhow::ensure!(
        temperature.iter().all(|t| t.is_finite() && *t > 0.0),
        "Emission temperatures must be positive and finite"
    );
    let integrated =
        inputs.spectral_integration_mode() == SpectralGridMode::AtmosphereIntegratedLineShape;
    let grid = if integrated {
        inputs.fine_spectral_grid()
    } else {
        inputs.spectral_grid()
    }
    .ok_or_else(|| anyhow!("Spectral grid must be set for line emission"))?;
    let wavenumbers = grid.central_wavenumber_cminv();
    let mut values = Array2::zeros((temperature.len(), wavenumbers.len()));
    let mut derivatives = d_areas.map(|_| Array2::zeros(values.raw_dim()));
    for (alt, &t) in temperature.iter().enumerate() {
        let mut row = values.row_mut(alt);
        for (i, &wavelength) in wavelengths.iter().enumerate() {
            let area = areas[[alt, i]];
            let center = 1.0e7 / wavelength;
            let width = OpticalLine::doppler_width_cminv(center, t, molecular_mass_amu);
            if let (Some(d_areas), Some(derivatives)) = (d_areas, derivatives.as_mut()) {
                let d_area = d_areas[[alt, i]];
                if area == 0.0 && d_area == 0.0 {
                    continue;
                }
                let line = AdjustedLineParameters {
                    line_center: center,
                    doppler_width: width,
                    line_intensity_re: area / (SQRT_PI * width),
                    ..Default::default()
                };
                let direction = LineShapeDirection {
                    doppler_width: width / (2.0 * t),
                    line_intensity_re: (d_area - area / (2.0 * t)) / (SQRT_PI * width),
                    ..Default::default()
                };
                assign_with_derivative(
                    wavenumbers.as_slice().unwrap(),
                    &line,
                    &direction,
                    LineShape::Gaussian,
                    row.as_slice_mut().unwrap(),
                    derivatives.row_mut(alt).as_slice_mut().unwrap(),
                );
            } else if area != 0.0 {
                assign_normalized_doppler_line_shape(
                    wavenumbers.as_slice().unwrap(),
                    center,
                    width,
                    area,
                    row.as_slice_mut().unwrap(),
                );
            }
        }
    }
    if grid.coordinate() == LineshapeCoordinate::WavelengthNm {
        for (i, &wavelength) in grid.central_wavelengths_nm().iter().enumerate() {
            let jacobian = 1.0e7 / wavelength.powi(2);
            values.column_mut(i).mapv_inplace(|v| v * jacobian);
            if let Some(d) = derivatives.as_mut() {
                d.column_mut(i).mapv_inplace(|v| v * jacobian);
            }
        }
    }
    if integrated {
        let mapping = inputs
            .spectral_mapping_matrix()
            .ok_or_else(|| anyhow!("Missing emission spectral mapping"))?;
        values = mapping.dot(values.view());
        derivatives = derivatives.map(|d| mapping.dot(d.view()));
    }
    Ok((values, derivatives))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atmosphere::types::SpectralGrid;
    use crate::emission::BandLineData;

    struct Inputs {
        altitude: Array1<f64>,
        temperature: Array1<f64>,
        grid: SpectralGrid,
    }

    impl Inputs {
        fn new() -> Self {
            Self {
                altitude: array![0.0, 1000.0, 2000.0],
                temperature: array![190.0, 240.0, 300.0],
                grid: SpectralGrid::from_wavenumbers_cminv(&Array1::linspace(
                    13099.75, 13100.75, 4001,
                ))
                .unwrap(),
            }
        }
    }

    impl StorageInputs for Inputs {
        fn num_stokes(&self) -> usize {
            1
        }
        fn spectral_integration_mode(&self) -> SpectralGridMode {
            SpectralGridMode::Monochromatic
        }
        fn num_singlescatter_moments(&self) -> usize {
            1
        }
        fn calculate_pressure_derivative(&self) -> bool {
            false
        }
        fn calculate_temperature_derivative(&self) -> bool {
            true
        }
        fn calculate_specific_humidity_derivative(&self) -> bool {
            false
        }
        fn altitude_m(&self) -> ArrayView1<'_, f64> {
            self.altitude.view()
        }
        fn pressure_pa(&self) -> Option<ArrayView1<'_, f64>> {
            None
        }
        fn temperature_k(&self) -> Option<ArrayView1<'_, f64>> {
            Some(self.temperature.view())
        }
        fn spectral_grid(&self) -> Option<&SpectralGrid> {
            Some(&self.grid)
        }
        fn fine_spectral_grid(&self) -> Option<&SpectralGrid> {
            Some(&self.grid)
        }
        fn air_numberdensity_dict(&self) -> HashMap<String, Array1<f64>> {
            HashMap::new()
        }
        fn dry_air_numberdensity_dict(&self) -> HashMap<String, Array1<f64>> {
            HashMap::new()
        }
    }

    /// A synthetic excitation model with no molecular level data or Boltzmann
    /// distribution. Its two weights vary with altitude and optionally kinetic T.
    struct SyntheticPopulation {
        lines: BandLineData,
        temperature_slope: f64,
        reject_derivatives: bool,
    }

    impl BandEmissionModel for SyntheticPopulation {
        fn line_data(&self) -> &BandLineData {
            &self.lines
        }

        fn line_weights(&self, inputs: &impl StorageInputs) -> Result<Array2<f64>> {
            let temperature = inputs.temperature_k().unwrap();
            let altitude = inputs.altitude_m();
            let mut weights = Array2::zeros((altitude.len(), 2));
            for i in 0..altitude.len() {
                let first = 0.25
                    + altitude[i] / 10_000.0
                    + self.temperature_slope * (temperature[i] - 230.0);
                weights[[i, 0]] = first;
                weights[[i, 1]] = 1.0 - first;
            }
            Ok(weights)
        }

        fn line_weights_with_temperature_derivative(
            &self,
            inputs: &impl StorageInputs,
        ) -> Result<(Array2<f64>, Array2<f64>)> {
            assert!(
                !self.reject_derivatives,
                "Value-only path requested derivatives"
            );
            let weights = self.line_weights(inputs)?;
            let mut derivatives = Array2::zeros(weights.raw_dim());
            derivatives.column_mut(0).fill(self.temperature_slope);
            derivatives.column_mut(1).fill(-self.temperature_slope);
            Ok((weights, derivatives))
        }
    }

    fn source(mass: f64, slope: f64) -> BandVolumeEmissionRate<SyntheticPopulation> {
        BandVolumeEmissionRate::new(
            array![0.0, 2000.0],
            array![0.0, 2.0],
            SyntheticPopulation {
                lines: BandLineData::new(
                    "synthetic",
                    array![1.0e7 / 13100.0, 1.0e7 / 13100.5],
                    mass,
                )
                .unwrap(),
                temperature_slope: slope,
                reject_derivatives: false,
            },
        )
        .unwrap()
    }

    fn integrate(values: ArrayView1<'_, f64>, grid: ArrayView1<'_, f64>) -> f64 {
        values
            .iter()
            .zip(values.iter().skip(1))
            .zip(grid.iter().zip(grid.iter().skip(1)))
            .map(|((&a, &b), (&x, &y))| (a + b) * (y - x) / 2.0)
            .sum()
    }

    #[test]
    fn arbitrary_populations_conserve_ver_and_have_correct_temperature_derivatives() {
        for slope in [0.0, 0.001] {
            let mut inputs = Inputs::new();
            let emission = source(16.0, slope);
            let ver = emission.interpolation(&inputs).dot(&emission.photon_ver);
            assert_eq!(ver, array![0.0, 1.0, 2.0]);
            let (value, derivative) = emission.spectrum(&inputs, true, Some(ver.view())).unwrap();
            let derivative = derivative.unwrap();
            let (unit_ver, _) = emission.spectrum(&inputs, false, None).unwrap();
            let grid = inputs.grid.central_wavenumber_cminv();
            for level in 0..3 {
                assert!((integrate(value.row(level), grid) - ver[level]).abs() < 1e-11);
                assert!(integrate(derivative.row(level), grid).abs() < 1e-11);
                assert!((integrate(unit_ver.row(level), grid) - 1.0).abs() < 1e-11);
            }

            let step = 0.001;
            inputs.temperature += step;
            let (above, _) = emission.spectrum(&inputs, false, Some(ver.view())).unwrap();
            inputs.temperature -= 2.0 * step;
            let (below, _) = emission.spectrum(&inputs, false, Some(ver.view())).unwrap();
            let numeric = (above - below) / (2.0 * step);
            let peak = derivative.iter().map(|v| v.abs()).fold(0.0, f64::max);
            assert!(
                peak > 0.0,
                "Fixed populations must still include Doppler d/dT"
            );
            for (analytic, numeric) in derivative.iter().zip(&numeric) {
                assert!((analytic - numeric).abs() < 1e-8 * peak);
            }
        }
    }

    #[test]
    fn molecular_mass_controls_width_without_changing_photon_area() {
        let inputs = Inputs::new();
        let (light, _) = source(16.0, 0.0).spectrum(&inputs, false, None).unwrap();
        let (heavy, _) = source(64.0, 0.0).spectrum(&inputs, false, None).unwrap();
        let grid = inputs.grid.central_wavenumber_cminv();
        for level in 0..3 {
            // Quadrupling mass halves the Doppler width and doubles each peak.
            let light_peak = light.row(level).iter().copied().fold(0.0, f64::max);
            let heavy_peak = heavy.row(level).iter().copied().fold(0.0, f64::max);
            assert!((heavy_peak / light_peak - 2.0).abs() < 1e-10);
            assert!((integrate(light.row(level), grid) - 1.0).abs() < 1e-11);
            assert!((integrate(heavy.row(level), grid) - 1.0).abs() < 1e-11);
        }
    }

    #[test]
    fn value_only_path_never_requests_population_derivatives() {
        let inputs = Inputs::new();
        let mut emission = source(16.0, 0.001);
        emission.model.reject_derivatives = true;
        let (_, derivative) = emission.spectrum(&inputs, false, None).unwrap();
        assert!(derivative.is_none());
    }

    #[test]
    fn invalid_spectroscopy_is_rejected() {
        for mass in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(BandLineData::new("test", array![760.0], mass).is_err());
        }
        for wavelengths in [array![], array![0.0], array![-1.0], array![f64::NAN]] {
            assert!(BandLineData::new("test", wavelengths, 16.0).is_err());
        }
    }
}
