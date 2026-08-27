//! Pure-Rust implementation of the MT_CKD 4.3 continuum.
//!
//! Pressure is in Pa, temperature in K, and volume mixing ratios are unitless.
//! The result is an extinction coefficient in m^-1 on the 1,991-point MT_CKD
//! grid from 0 through 19,900 cm^-1. All five atmospheric derivatives are
//! propagated analytically in the same evaluation when requested.

use std::fmt;
use std::fs;
use std::io;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::path::Path;

/// Number of useful MT_CKD samples from 0 through 19,900 cm^-1.
pub const SPECTRAL_POINTS: usize = 1991;
/// Number of values returned by the `mt_ckd` interface.
pub const OUTPUT_POINTS: usize = SPECTRAL_POINTS;
/// Wavenumber represented by output index zero.
pub const OUTPUT_WAVENUMBER_START: f64 = 0.0;
/// Fixed-grid spacing in cm^-1.
pub const WAVENUMBER_SPACING: f64 = 10.0;

const N_INPUTS: usize = 5;
const N_FIELDS: usize = 22;
const HEADER_BYTES: usize = 16;
const RADCN2: f64 = 1.438_775_2;
const ALOSMT: f64 = 2.686_777_5e19;
const XLOSMT: f64 = 2.686_75e19;
// MT_CKD computes optical depth from a column amount. A one-metre column makes
// that optical depth numerically equal to the extinction coefficient in m^-1.
const REFERENCE_PATH_LENGTH_CM: f64 = 100.0;
const N2_FUNDAMENTAL_POINTS: usize = 228;
const N2_FUNDAMENTAL_V1: f64 = 1_997.784_896;
const N2_FUNDAMENTAL_DV: f64 = 3.981_461_525;
const TABLE_MAGIC: &[u8; 8] = b"MTCKD43\0";

#[derive(Clone, Copy, Debug)]
#[repr(usize)]
enum Input {
    Pressure = 0,
    Temperature = 1,
    H2O = 2,
    CO2 = 3,
    O3 = 4,
}

#[derive(Clone, Copy, Debug, Default)]
struct Dual {
    value: f64,
    derivative: [f64; N_INPUTS],
}

trait Scalar:
    Copy
    + Add<Output = Self>
    + Add<f64, Output = Self>
    + Sub<Output = Self>
    + Sub<f64, Output = Self>
    + Mul<Output = Self>
    + Mul<f64, Output = Self>
    + Div<Output = Self>
    + Div<f64, Output = Self>
    + Neg<Output = Self>
{
    fn constant(value: f64) -> Self;
    fn value(self) -> f64;
    fn exp(self) -> Self;
    fn powf(self, exponent: f64) -> Self;
}

impl Scalar for f64 {
    fn constant(value: f64) -> Self {
        value
    }

    fn value(self) -> f64 {
        self
    }

    fn exp(self) -> Self {
        f64::exp(self)
    }

    fn powf(self, exponent: f64) -> Self {
        f64::powf(self, exponent)
    }
}

impl Dual {
    fn constant(value: f64) -> Self {
        Self {
            value,
            derivative: [0.0; N_INPUTS],
        }
    }

    fn variable(value: f64, input: Input) -> Self {
        let mut derivative = [0.0; N_INPUTS];
        derivative[input as usize] = 1.0;
        Self { value, derivative }
    }

    fn exp(self) -> Self {
        let value = self.value.exp();
        Self {
            value,
            derivative: self.derivative.map(|d| d * value),
        }
    }

    fn powf(self, exponent: f64) -> Self {
        let value = self.value.powf(exponent);
        let scale = exponent * self.value.powf(exponent - 1.0);
        Self {
            value,
            derivative: self.derivative.map(|d| d * scale),
        }
    }
}

impl Scalar for Dual {
    fn constant(value: f64) -> Self {
        Self::constant(value)
    }

    fn value(self) -> f64 {
        self.value
    }

    fn exp(self) -> Self {
        self.exp()
    }

    fn powf(self, exponent: f64) -> Self {
        self.powf(exponent)
    }
}

impl Add for Dual {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value + rhs.value,
            derivative: std::array::from_fn(|i| self.derivative[i] + rhs.derivative[i]),
        }
    }
}

impl Add<f64> for Dual {
    type Output = Self;

    fn add(self, rhs: f64) -> Self::Output {
        self + Self::constant(rhs)
    }
}

impl Sub for Dual {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value - rhs.value,
            derivative: std::array::from_fn(|i| self.derivative[i] - rhs.derivative[i]),
        }
    }
}

impl Sub<f64> for Dual {
    type Output = Self;

    fn sub(self, rhs: f64) -> Self::Output {
        self - Self::constant(rhs)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Mul for Dual {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value * rhs.value,
            derivative: std::array::from_fn(|i| {
                self.derivative[i] * rhs.value + self.value * rhs.derivative[i]
            }),
        }
    }
}

impl Mul<f64> for Dual {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            value: self.value * rhs,
            derivative: self.derivative.map(|d| d * rhs),
        }
    }
}

impl Div for Dual {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let denominator = rhs.value * rhs.value;
        Self {
            value: self.value / rhs.value,
            derivative: std::array::from_fn(|i| {
                (self.derivative[i] * rhs.value - self.value * rhs.derivative[i]) / denominator
            }),
        }
    }
}

impl Div<f64> for Dual {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        self * (1.0 / rhs)
    }
}

impl Neg for Dual {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            value: -self.value,
            derivative: self.derivative.map(|d| -d),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct AtmosphericState {
    pub pressure_pa: f64,
    pub temperature_k: f64,
    pub h2o_vmr: f64,
    pub co2_vmr: f64,
    pub o3_vmr: f64,
}

/// MT_CKD extinction coefficient and analytic Jacobians on the fixed grid.
#[derive(Debug, Clone)]
pub struct LinearizedSpectrum {
    pub extinction: Vec<f64>,
    pub d_pressure_pa: Vec<f64>,
    pub d_temperature_k: Vec<f64>,
    pub d_h2o_vmr: Vec<f64>,
    pub d_co2_vmr: Vec<f64>,
    pub d_o3_vmr: Vec<f64>,
}

#[derive(Clone, Copy)]
#[repr(usize)]
enum Field {
    H2OSelf = 0,
    H2OSelfTemperatureExponent = 1,
    H2OForeign = 2,
    CO2 = 3,
    CO2TemperatureExponent = 4,
    O3Constant = 5,
    O3Linear = 6,
    O3Quadratic = 7,
    O2Fundamental = 8,
    O2FundamentalEnergy = 9,
    O2Infrared1 = 10,
    O2Infrared2 = 11,
    O2Infrared3 = 12,
    O2Visible = 13,
    N2Rotation296 = 14,
    N2RotationScale296 = 15,
    N2Rotation220 = 16,
    N2RotationScale220 = 17,
    N2Fundamental272 = 18,
    N2Fundamental228 = 19,
    N2H2OEfficiency = 20,
    N2Overtone = 21,
}

/// Errors raised while loading the separately distributed MT_CKD data table.
#[derive(Debug)]
pub enum MtCkdDataError {
    Io(io::Error),
    InvalidFormat(String),
}

impl fmt::Display for MtCkdDataError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(error) => write!(formatter, "failed to read MT_CKD data: {error}"),
            Self::InvalidFormat(message) => write!(formatter, "invalid MT_CKD data: {message}"),
        }
    }
}

impl std::error::Error for MtCkdDataError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            Self::InvalidFormat(_) => None,
        }
    }
}

impl From<io::Error> for MtCkdDataError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

/// An MT_CKD 4.3 evaluator loaded from the separately licensed coefficient data.
pub struct MtCkd {
    coefficients: Box<[f64]>,
}

impl MtCkd {
    /// Load and validate an MT_CKD 4.3 coefficient table from memory.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, MtCkdDataError> {
        if bytes.len() < HEADER_BYTES {
            return Err(MtCkdDataError::InvalidFormat(format!(
                "header is truncated: expected at least {HEADER_BYTES} bytes, found {}",
                bytes.len()
            )));
        }
        if &bytes[..TABLE_MAGIC.len()] != TABLE_MAGIC {
            return Err(MtCkdDataError::InvalidFormat(
                "unexpected magic or model version".to_owned(),
            ));
        }

        let spectral_points = u32::from_le_bytes(bytes[8..12].try_into().unwrap()) as usize;
        let fields = u32::from_le_bytes(bytes[12..16].try_into().unwrap()) as usize;
        if spectral_points != SPECTRAL_POINTS || fields != N_FIELDS {
            return Err(MtCkdDataError::InvalidFormat(format!(
                "expected {SPECTRAL_POINTS} spectral points and {N_FIELDS} fields, found {spectral_points} and {fields}"
            )));
        }

        let expected_length = HEADER_BYTES + SPECTRAL_POINTS * N_FIELDS * size_of::<f64>();
        if bytes.len() != expected_length {
            return Err(MtCkdDataError::InvalidFormat(format!(
                "expected {expected_length} bytes, found {}",
                bytes.len()
            )));
        }

        let coefficients: Box<[f64]> = bytes[HEADER_BYTES..]
            .chunks_exact(size_of::<f64>())
            .map(|value| f64::from_le_bytes(value.try_into().unwrap()))
            .collect();
        if let Some(index) = coefficients.iter().position(|value| !value.is_finite()) {
            return Err(MtCkdDataError::InvalidFormat(format!(
                "coefficient {index} is not finite"
            )));
        }

        Ok(Self { coefficients })
    }

    /// Load and validate an MT_CKD 4.3 coefficient table from a file.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self, MtCkdDataError> {
        Self::from_bytes(&fs::read(path)?)
    }

    fn coefficient(&self, field: Field, spectral_index: usize) -> f64 {
        self.coefficients[field as usize * SPECTRAL_POINTS + spectral_index]
    }
}

fn radiation_term<S: Scalar>(wavenumber: f64, temperature: S) -> S {
    let ratio = S::constant(wavenumber * RADCN2) / temperature;
    if ratio.value() <= 0.01 {
        ratio * (0.5 * wavenumber)
    } else if ratio.value() <= 10.0 {
        let exponential = (-ratio).exp();
        S::constant(wavenumber) * (S::constant(1.0) - exponential)
            / (S::constant(1.0) + exponential)
    } else {
        S::constant(wavenumber)
    }
}

fn interpolate_logarithmically<S: Scalar>(
    at_first_temperature: f64,
    at_second_temperature: f64,
    fraction: S,
) -> S {
    if at_first_temperature > 0.0 && at_second_temperature > 0.0 {
        S::constant(at_first_temperature)
            * (fraction * (at_second_temperature / at_first_temperature).ln()).exp()
    } else {
        S::constant(at_first_temperature)
            + fraction * (at_second_temperature - at_first_temperature)
    }
}

fn n2_fundamental_source<S: Scalar>(model: &MtCkd, position: isize, temperature: S) -> [S; 3] {
    // n2_ver_1 pads two zero samples before and after its stored grid.
    let raw_index = position - 3;
    if raw_index < 0 || raw_index >= N2_FUNDAMENTAL_POINTS as isize {
        return [S::constant(0.0); 3];
    }
    let raw_index = raw_index as usize;
    let at_272 = model.coefficient(Field::N2Fundamental272, raw_index);
    let at_228 = model.coefficient(Field::N2Fundamental228, raw_index);
    let reciprocal_temperature_fraction =
        (S::constant(1.0) / temperature - 1.0 / 272.0) / ((1.0 / 228.0) - (1.0 / 272.0));
    let linear_temperature_fraction = (temperature - 272.0) / (228.0 - 272.0);
    let spectral_coefficient = if at_272 > 0.0 && at_228 > 0.0 {
        interpolate_logarithmically(at_272, at_228, reciprocal_temperature_fraction)
    } else {
        S::constant(at_272) + linear_temperature_fraction * (at_228 - at_272)
    };
    let source_wavenumber = N2_FUNDAMENTAL_V1 + raw_index as f64 * N2_FUNDAMENTAL_DV;
    let spectral_coefficient = spectral_coefficient / source_wavenumber;
    let o2_efficiency = S::constant(1.294) - temperature * (0.4545 / 296.0);
    let h2o_efficiency = model.coefficient(Field::N2H2OEfficiency, raw_index);
    [
        spectral_coefficient,
        spectral_coefficient * o2_efficiency,
        spectral_coefficient * (9.0 / 7.0 * h2o_efficiency),
    ]
}

fn n2_fundamental_coefficients<S: Scalar>(
    model: &MtCkd,
    wavenumber: f64,
    temperature: S,
) -> [S; 3] {
    if !(2000.0..=2900.0).contains(&wavenumber) {
        return [S::constant(0.0); 3];
    }

    let padded_v1 = N2_FUNDAMENTAL_V1 - 2.0 * N2_FUNDAMENTAL_DV;
    let position = ((wavenumber - padded_v1) / N2_FUNDAMENTAL_DV + 1.001).trunc() as isize;
    let position_wavenumber = padded_v1 + (position - 1) as f64 * N2_FUNDAMENTAL_DV;
    let p = (wavenumber - position_wavenumber) / N2_FUNDAMENTAL_DV;
    let c = (3.0 - 2.0 * p) * p * p;
    let b = 0.5 * p * (1.0 - p);
    let b1 = b * (1.0 - p);
    let b2 = b * p;
    let weights = [-b1, 1.0 - c + b2, c + b1, -b2];
    let mut result = [S::constant(0.0); 3];
    for (offset, weight) in weights.into_iter().enumerate() {
        let source = n2_fundamental_source(model, position - 1 + offset as isize, temperature);
        for component in 0..3 {
            result[component] = result[component] + source[component] * weight;
        }
    }
    result
}

fn evaluate<S: Scalar, const INCLUDE_RAYLEIGH: bool>(
    model: &MtCkd,
    pressure_pa: S,
    temperature: S,
    h2o_vmr: S,
    co2_vmr: S,
    o3_vmr: S,
) -> Vec<S> {
    let pressure_mb = pressure_pa / 100.0;
    let initial_total = pressure_mb / 1013.0
        * (S::constant(273.0) / temperature)
        * (ALOSMT * REFERENCE_PATH_LENGTH_CM);
    let dry_amount = initial_total * (S::constant(1.0) - h2o_vmr);
    // Preserve MT_CKD's special pure-water convention. The selected branch is
    // differentiated analytically; as in the upstream piecewise model, the
    // derivative is undefined exactly at either threshold.
    let h2o_amount = if (h2o_vmr.value() - 1.0).abs() < 1.0e-5 {
        initial_total
    } else {
        h2o_vmr * dry_amount
    };
    let co2_amount = co2_vmr * dry_amount;
    let o3_amount = o3_vmr * dry_amount;
    let o2_amount = dry_amount * 0.21;
    let broadening_amount = dry_amount * (0.78 + 0.009);
    let total_amount = broadening_amount + h2o_amount + co2_amount + o3_amount + o2_amount;
    let x_h2o = h2o_amount / total_amount;
    let x_o2 = o2_amount / total_amount;
    let x_n2 = S::constant(1.0) - x_h2o - x_o2;
    let n2_amount = x_n2 * total_amount;
    let rho = pressure_mb / 1013.0 * (S::constant(296.0) / temperature);
    let amagat = pressure_mb / 1013.0 * (S::constant(273.0) / temperature);

    let mut extinction = Vec::with_capacity(OUTPUT_POINTS);

    for spectral_index in 0..SPECTRAL_POINTS {
        let wavenumber = spectral_index as f64 * WAVENUMBER_SPACING;
        let radiation = radiation_term(wavenumber, temperature);

        let self_coefficient = S::constant(model.coefficient(Field::H2OSelf, spectral_index))
            * (S::constant(296.0) / temperature)
                .powf(model.coefficient(Field::H2OSelfTemperatureExponent, spectral_index))
            * x_h2o
            * rho;
        let foreign_coefficient = S::constant(model.coefficient(Field::H2OForeign, spectral_index))
            * (S::constant(1.0) - x_h2o)
            * rho;
        let h2o = h2o_amount * (self_coefficient + foreign_coefficient);

        let co2 = co2_amount
            * rho
            * 1.0e-20
            * model.coefficient(Field::CO2, spectral_index)
            * (temperature / 246.0)
                .powf(model.coefficient(Field::CO2TemperatureExponent, spectral_index));

        let delta_temperature = temperature - 273.15;
        let o3_coefficient = S::constant(model.coefficient(Field::O3Constant, spectral_index))
            + delta_temperature * model.coefficient(Field::O3Linear, spectral_index)
            + delta_temperature
                * delta_temperature
                * model.coefficient(Field::O3Quadratic, spectral_index);
        let o3 = o3_amount * 1.0e-20 * o3_coefficient;

        let o2_fundamental = o2_amount
            * 1.0e-20
            * amagat
            * model.coefficient(Field::O2Fundamental, spectral_index)
            * (S::constant(model.coefficient(Field::O2FundamentalEnergy, spectral_index))
                * (S::constant(1.0 / 296.0) - S::constant(1.0) / temperature))
                .exp();
        let o2_infrared1 = o2_amount / XLOSMT
            * amagat
            * (x_o2 / 0.446 + x_n2 * (0.3 / 0.446) + x_h2o)
            * model.coefficient(Field::O2Infrared1, spectral_index);
        let o2_infrared2 = o2_amount / total_amount
            * (1.0 / 0.209)
            * o2_amount
            * 1.0e-20
            * rho
            * model.coefficient(Field::O2Infrared2, spectral_index);
        let o2_infrared3 =
            o2_amount / XLOSMT * amagat * model.coefficient(Field::O2Infrared3, spectral_index);
        let o2_visible = o2_amount / total_amount
            * o2_amount
            * 1.0e-20
            * amagat
            * model.coefficient(Field::O2Visible, spectral_index);

        let rotation_fraction = (temperature - 296.0) / (220.0 - 296.0);
        let n2_rotation_coefficient = interpolate_logarithmically(
            model.coefficient(Field::N2Rotation296, spectral_index),
            model.coefficient(Field::N2Rotation220, spectral_index),
            rotation_fraction,
        );
        let n2_rotation_scale = interpolate_logarithmically(
            model.coefficient(Field::N2RotationScale296, spectral_index),
            model.coefficient(Field::N2RotationScale220, spectral_index),
            rotation_fraction,
        );
        let n2_o2_efficiency = (n2_rotation_scale - 1.0) * (0.79 / 0.21);
        let n2_rotation = n2_amount / XLOSMT
            * amagat
            * n2_rotation_coefficient
            * (x_n2 + n2_o2_efficiency * x_o2 + x_h2o);

        let n2_fundamental_coefficients =
            n2_fundamental_coefficients(model, wavenumber, temperature);
        let n2_fundamental = n2_amount / XLOSMT
            * amagat
            * (x_n2 * n2_fundamental_coefficients[0]
                + x_o2 * n2_fundamental_coefficients[1]
                + x_h2o * n2_fundamental_coefficients[2]);
        let n2_overtone = n2_amount / XLOSMT
            * amagat
            * (x_n2 + x_o2 + x_h2o)
            * model.coefficient(Field::N2Overtone, spectral_index);

        let absorption = (h2o
            + co2
            + o3
            + o2_fundamental
            + o2_infrared1
            + o2_infrared2
            + o2_infrared3
            + o2_visible
            + n2_rotation
            + n2_fundamental
            + n2_overtone)
            * radiation;

        // MT_CKD's standalone wrapper also enables its historical Rayleigh
        // term. With the radiation term restored this simplifies exactly to
        // the expression below.
        let scaled_wavenumber = wavenumber / 1.0e4;
        let rayleigh = total_amount
            * (1.0e-20 / (2.686_75e-1 * 1.0e5))
            * (scaled_wavenumber.powi(4) / (9.380_76e2 - 10.8426 * scaled_wavenumber.powi(2)));

        let result = if INCLUDE_RAYLEIGH {
            absorption + rayleigh
        } else {
            absorption
        };
        extinction.push(result);
    }

    extinction
}

fn calculate_linearized_impl<const INCLUDE_RAYLEIGH: bool>(
    model: &MtCkd,
    state: AtmosphericState,
) -> LinearizedSpectrum {
    let spectrum = evaluate::<Dual, INCLUDE_RAYLEIGH>(
        model,
        Dual::variable(state.pressure_pa, Input::Pressure),
        Dual::variable(state.temperature_k, Input::Temperature),
        Dual::variable(state.h2o_vmr, Input::H2O),
        Dual::variable(state.co2_vmr, Input::CO2),
        Dual::variable(state.o3_vmr, Input::O3),
    );

    let mut extinction = Vec::with_capacity(OUTPUT_POINTS);
    let mut derivatives =
        std::array::from_fn::<_, N_INPUTS, _>(|_| Vec::with_capacity(OUTPUT_POINTS));
    for result in spectrum {
        extinction.push(result.value);
        for (input, output) in derivatives.iter_mut().enumerate() {
            output.push(result.derivative[input]);
        }
    }

    LinearizedSpectrum {
        extinction,
        d_pressure_pa: std::mem::take(&mut derivatives[Input::Pressure as usize]),
        d_temperature_k: std::mem::take(&mut derivatives[Input::Temperature as usize]),
        d_h2o_vmr: std::mem::take(&mut derivatives[Input::H2O as usize]),
        d_co2_vmr: std::mem::take(&mut derivatives[Input::CO2 as usize]),
        d_o3_vmr: std::mem::take(&mut derivatives[Input::O3 as usize]),
    }
}

/// Evaluate the complete MT_CKD 4.3 standalone spectrum, including its
/// historical Rayleigh-scattering term, and all atmospheric Jacobians.
pub fn calculate_linearized(model: &MtCkd, state: AtmosphericState) -> LinearizedSpectrum {
    calculate_linearized_impl::<true>(model, state)
}

/// Evaluate MT_CKD 4.3 continuum absorption and all atmospheric Jacobians,
/// excluding the historical Rayleigh-scattering term.
pub fn calculate_absorption_linearized(
    model: &MtCkd,
    state: AtmosphericState,
) -> LinearizedSpectrum {
    calculate_linearized_impl::<false>(model, state)
}

/// Evaluate the complete MT_CKD 4.3 standalone spectrum, including its
/// historical Rayleigh-scattering term, without derivative propagation.
pub fn calculate(model: &MtCkd, state: AtmosphericState) -> Vec<f64> {
    evaluate::<f64, true>(
        model,
        state.pressure_pa,
        state.temperature_k,
        state.h2o_vmr,
        state.co2_vmr,
        state.o3_vmr,
    )
}

/// Evaluate only MT_CKD 4.3 continuum absorption, excluding the historical
/// Rayleigh-scattering term, without derivative propagation.
pub fn calculate_absorption(model: &MtCkd, state: AtmosphericState) -> Vec<f64> {
    evaluate::<f64, false>(
        model,
        state.pressure_pa,
        state.temperature_k,
        state.h2o_vmr,
        state.co2_vmr,
        state.o3_vmr,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn state() -> AtmosphericState {
        AtmosphericState {
            pressure_pa: 101_325.0,
            temperature_k: 288.15,
            h2o_vmr: 0.01,
            co2_vmr: 400.0e-6,
            o3_vmr: 2.0e-6,
        }
    }

    fn zero_coefficient_model() -> MtCkd {
        let mut bytes =
            Vec::with_capacity(HEADER_BYTES + SPECTRAL_POINTS * N_FIELDS * size_of::<f64>());
        bytes.extend_from_slice(TABLE_MAGIC);
        bytes.extend_from_slice(&(SPECTRAL_POINTS as u32).to_le_bytes());
        bytes.extend_from_slice(&(N_FIELDS as u32).to_le_bytes());
        bytes.resize(
            HEADER_BYTES + SPECTRAL_POINTS * N_FIELDS * size_of::<f64>(),
            0,
        );
        MtCkd::from_bytes(&bytes).unwrap()
    }

    #[test]
    fn coefficient_table_parser_validates_the_data_contract() {
        let model = zero_coefficient_model();
        assert_eq!(model.coefficients.len(), N_FIELDS * SPECTRAL_POINTS);

        let error = match MtCkd::from_bytes(b"not an MT_CKD table") {
            Ok(_) => panic!("invalid MT_CKD data was accepted"),
            Err(error) => error,
        };
        assert!(error.to_string().contains("unexpected magic"));
    }

    #[test]
    fn analytic_derivatives_match_central_differences() {
        let model = zero_coefficient_model();
        let reference = calculate_linearized(&model, state());
        let perturbations = [1.0, 1.0e-3, 1.0e-7, 1.0e-8, 1.0e-8];
        let analytic = [
            &reference.d_pressure_pa,
            &reference.d_temperature_k,
            &reference.d_h2o_vmr,
            &reference.d_co2_vmr,
            &reference.d_o3_vmr,
        ];

        for input in 0..N_INPUTS {
            let mut above = state();
            let mut below = state();
            let delta = perturbations[input];
            match input {
                0 => {
                    above.pressure_pa += delta;
                    below.pressure_pa -= delta;
                }
                1 => {
                    above.temperature_k += delta;
                    below.temperature_k -= delta;
                }
                2 => {
                    above.h2o_vmr += delta;
                    below.h2o_vmr -= delta;
                }
                3 => {
                    above.co2_vmr += delta;
                    below.co2_vmr -= delta;
                }
                4 => {
                    above.o3_vmr += delta;
                    below.o3_vmr -= delta;
                }
                _ => unreachable!(),
            }
            let above = calculate(&model, above);
            let below = calculate(&model, below);
            for spectral_index in 1..SPECTRAL_POINTS {
                let numeric = (above[spectral_index] - below[spectral_index]) / (2.0 * delta);
                let scale = numeric.abs().max(analytic[input][spectral_index].abs());
                if scale > 1.0e-18 {
                    assert!(
                        (numeric - analytic[input][spectral_index]).abs() <= 2.0e-5 * scale,
                        "input={input}, index={spectral_index}, numeric={numeric:e}, analytic={:e}",
                        analytic[input][spectral_index]
                    );
                }
            }
        }
    }

    #[test]
    fn value_only_evaluation_matches_linearized_values() {
        let model = zero_coefficient_model();
        assert_eq!(
            calculate(&model, state()),
            calculate_linearized(&model, state()).extinction
        );
        assert_eq!(
            calculate_absorption(&model, state()),
            calculate_absorption_linearized(&model, state()).extinction
        );
    }

    #[test]
    fn pure_water_branch_is_finite_and_linearized() {
        let model = zero_coefficient_model();
        let mut pure_water = state();
        pure_water.h2o_vmr = 1.0;

        let reference = calculate_linearized(&model, pure_water);
        assert!(reference.extinction.iter().all(|value| value.is_finite()));
        assert!(reference.d_h2o_vmr.iter().all(|value| value.is_finite()));

        let delta = 1.0e-7;
        let mut above = pure_water;
        let mut below = pure_water;
        above.h2o_vmr += delta;
        below.h2o_vmr -= delta;
        let above = calculate(&model, above);
        let below = calculate(&model, below);
        for spectral_index in 1..SPECTRAL_POINTS {
            let numeric = (above[spectral_index] - below[spectral_index]) / (2.0 * delta);
            let analytic = reference.d_h2o_vmr[spectral_index];
            let scale = numeric.abs().max(analytic.abs());
            if scale > 1.0e-18 {
                assert!(
                    (numeric - analytic).abs() <= 2.0e-7 * scale,
                    "index={spectral_index}, numeric={numeric:e}, analytic={analytic:e}"
                );
            }
        }
    }
}
