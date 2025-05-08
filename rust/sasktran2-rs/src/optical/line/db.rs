use crate::prelude::*;

use std::f64::consts::FRAC_1_SQRT_2;

const C2: f64 = 1.4387769;
const SPEED_OF_LIGHT: f64 = 2.99792458e10;
const NA: f64 = 6.02214179e23;
const K_B: f64 = 1.38064852e-16;
const SQRT_PI: f64 = 1.772_453_850_905_516;

pub struct OpticalLine {
    pub line_center: f64,
    pub line_intensity: f64,
    pub lower_energy: f64,
    pub gamma_air: f64,
    pub gamma_self: f64,
    pub delta_air: f64,
    pub n_air: f64,
    pub mol_id: i32,
    pub iso_id: i32,
    pub y_coupling: Vec<f64>,
    pub g_coupling: Vec<f64>,
    pub coupling_temperature: Vec<f64>,
}

pub struct AdjustedLineParameters {
    pub line_center: f64,
    pub line_intensity_re: f64,
    pub line_intensity_im: f64,
    pub doppler_width: f64,
    pub y: f64,
}

pub struct OpticalLineDB {
    pub lines: Vec<OpticalLine>,
}

impl OpticalLine {
    #[inline(always)]
    pub fn adjusted_parameters(
        &self,
        temperature: f64,
        pressure: f64,
        pself: f64,
        partition_factor: f64,
        mol_mass: f64,
    ) -> Result<AdjustedLineParameters> {
        // in atm
        let pressure = pressure / 101325.0;
        let pself = pself / 101325.0;

        let le = self.lower_energy;
        let lc = self.line_center;

        let doppler_width =
            lc / SPEED_OF_LIGHT * (NA * K_B * temperature / mol_mass).sqrt() / FRAC_1_SQRT_2;

        let numerator = (-C2 * le / temperature).exp() * (1.0 - (-C2 * lc / temperature).exp());

        let denominator = (1.0 - (-C2 * lc / 296.0).exp()) * (-C2 * le / 296.0).exp();

        let mut line_intensity_re =
            self.line_intensity * numerator / denominator / partition_factor * 1.0
                / (SQRT_PI * doppler_width);

        line_intensity_re /= 1e4; // to m^2

        let da = self.delta_air;
        let lc = self.line_center;

        let line_center = lc + da * pressure;

        let na = self.n_air;
        let g_air = self.gamma_air;
        let g_self = self.gamma_self;
        let gamma_val =
            (296.0 / temperature).powf(na) * (g_air * (pressure - pself) + g_self * pself);

        let y = gamma_val / doppler_width;

        Ok(AdjustedLineParameters {
            line_center,
            line_intensity_re,
            line_intensity_im: 0.0,
            doppler_width,
            y,
        })
    }
}

impl OpticalLineDB {
    pub fn between_slice(&self, min: f64, max: f64) -> &[OpticalLine] {
        let mut start = 0;
        let mut end = self.lines.len();

        for (i, line) in self.lines.iter().enumerate() {
            if line.line_center > min {
                start = i;
                break;
            }
        }

        for (i, line) in self.lines.iter().enumerate() {
            if line.line_center > max {
                end = i;
                break;
            }
        }

        &self.lines[start..end]
    }
}
