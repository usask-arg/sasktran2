use crate::optical::storage::OpticalQuantities;
use anyhow::*;
use gauss_quad::GaussLegendre;
use ndarray::Array1;

use super::StorageInputs;

/// Object allowing the inputs for the storage to be set manually
/// This is typically used by functions like optical_property.cross_section that wants to piggy back on
/// optical_property.atmosphere_quantities(atmo)
#[derive(Default)]
#[allow(dead_code)]
pub struct ManualStorageInputs {
    num_stokes: Option<usize>,
    num_singlescatter_moments: Option<usize>,
    calculate_pressure_derivative: Option<bool>,
    calculate_temperature_derivative: Option<bool>,
    calculate_specific_humidity_derivative: Option<bool>,
    altitude_m: Option<ndarray::Array1<f64>>,
    pressure_pa: Option<ndarray::Array1<f64>>,
    temperature_k: Option<ndarray::Array1<f64>>,
    wavelengths_nm: Option<ndarray::Array1<f64>>,
    wavenumbers_cminv: Option<ndarray::Array1<f64>>,
    air_numberdensity_dict: std::collections::HashMap<String, ndarray::Array1<f64>>,
    dry_air_numberdensity_dict: std::collections::HashMap<String, ndarray::Array1<f64>>,
}

impl ManualStorageInputs {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_num_stokes(mut self, num_stokes: usize) -> Self {
        self.num_stokes = Some(num_stokes);
        self
    }

    pub fn with_wavelengths_nm(mut self, wavelengths_nm: ndarray::Array1<f64>) -> Self {
        self.wavenumbers_cminv = Some(wavelengths_nm.clone().mapv(|x| 1.0e7 / x));

        self.wavelengths_nm = Some(wavelengths_nm);
        self
    }

    pub fn with_singlescatter_moments(mut self, num_singlescatter_moments: usize) -> Self {
        self.num_singlescatter_moments = Some(num_singlescatter_moments);
        self
    }

    pub fn with_altitude_m(mut self, altitude_m: ndarray::Array1<f64>) -> Self {
        self.altitude_m = Some(altitude_m);
        self
    }
}

impl StorageInputs for ManualStorageInputs {
    fn num_stokes(&self) -> usize {
        self.num_stokes.expect("num_stokes not set")
    }

    fn num_singlescatter_moments(&self) -> usize {
        self.num_singlescatter_moments
            .expect("num_singlescatter_moments not set")
    }

    fn calculate_pressure_derivative(&self) -> bool {
        todo!()
    }

    fn calculate_temperature_derivative(&self) -> bool {
        todo!()
    }

    fn calculate_specific_humidity_derivative(&self) -> bool {
        todo!()
    }

    fn altitude_m(&self) -> ndarray::ArrayView1<'_, f64> {
        self.altitude_m.as_ref().expect("altitude_m not set").view()
    }

    fn pressure_pa(&self) -> Option<ndarray::ArrayView1<'_, f64>> {
        todo!()
    }

    fn temperature_k(&self) -> Option<ndarray::ArrayView1<'_, f64>> {
        todo!()
    }

    fn wavelengths_nm(&self) -> Option<ndarray::ArrayView1<'_, f64>> {
        self.wavelengths_nm.as_ref().map(|x| x.view())
    }

    fn wavenumbers_cminv(&self) -> Option<ndarray::ArrayView1<'_, f64>> {
        self.wavenumbers_cminv.as_ref().map(|x| x.view())
    }

    fn air_numberdensity_dict(&self) -> std::collections::HashMap<String, ndarray::Array1<f64>> {
        todo!()
    }

    fn dry_air_numberdensity_dict(
        &self,
    ) -> std::collections::HashMap<String, ndarray::Array1<f64>> {
        todo!()
    }

    fn wavenumbers_cminv_left(&self) -> Option<ndarray::ArrayView1<'_, f64>> {
        todo!()
    }

    fn wavenumbers_cminv_right(&self) -> Option<ndarray::ArrayView1<'_, f64>> {
        todo!()
    }
}

pub struct StorageInputsQuadrature<'a> {
    base_input: &'a dyn StorageInputs,
    hires_wvnum: Array1<f64>,
    quadrature_x: Vec<f64>,
    quadrature_w: Vec<f64>,
}

impl StorageInputsQuadrature<'_> {
    pub fn from_base_input(base_input: &dyn StorageInputs) -> StorageInputsQuadrature<'_> {
        StorageInputsQuadrature {
            base_input,
            hires_wvnum: Array1::<f64>::zeros(0),
            quadrature_x: vec![0.0],
            quadrature_w: vec![1.0],
        }
    }

    pub fn with_quadrature_order(mut self, order: usize) -> Result<Self> {
        let quad = GaussLegendre::new(order)?;

        let mut nodes = vec![];
        let mut weights = vec![];

        for (node, w) in quad.iter().rev() {
            nodes.push((*node + 1.0) / 2.0);
            weights.push(*w / 2.0);
        }
        self.quadrature_w = weights;
        self.quadrature_x = nodes;

        let num_base = self.base_input.wavenumbers_cminv().unwrap().len();

        // Construct a new wavenumber grid
        let mut wavenumbers_cminv = Array1::<f64>::zeros(num_base * order);
        let wavenum_left = self.base_input.wavenumbers_cminv_left().unwrap();
        let wavenum_right = self.base_input.wavenumbers_cminv_right().unwrap();

        for i in 0..num_base {
            let left = wavenum_left[i];
            let right = wavenum_right[i];

            let delta = right - left;

            for j in 0..order {
                wavenumbers_cminv[order * i + j] = left + self.quadrature_x[j] * delta;
            }
        }
        self.hires_wvnum = wavenumbers_cminv;

        Ok(self)
    }

    pub fn reduce(&self, oq: &OpticalQuantities) -> Result<OpticalQuantities> {
        // Apply quadrature reduction to OpticalQuantities
        let mut oq_reduced = OpticalQuantities::default();

        let num_result = self
            .base_input
            .wavenumbers_cminv_left()
            .ok_or_else(|| anyhow!("wavenumbers_cminv_left not set"))?
            .len();

        let num_geo = oq.cross_section.shape()[0];

        oq_reduced.resize(num_geo, num_result);

        for i in 0..num_result {
            for j in 0..num_geo {
                oq_reduced.cross_section[[j, i]] = self
                    .quadrature_w
                    .iter()
                    .enumerate()
                    .map(|(k, w)| oq.cross_section[[j, self.quadrature_x.len() * i + k]] * w)
                    .sum();
                oq_reduced.ssa[[j, i]] = self
                    .quadrature_w
                    .iter()
                    .enumerate()
                    .map(|(k, w)| oq.ssa[[j, self.quadrature_x.len() * i + k]] * w)
                    .sum();
            }
        }

        Ok(oq_reduced)
    }
}

impl StorageInputs for StorageInputsQuadrature<'_> {
    fn num_stokes(&self) -> usize {
        self.base_input.num_stokes()
    }

    fn num_singlescatter_moments(&self) -> usize {
        self.base_input.num_singlescatter_moments()
    }

    fn calculate_pressure_derivative(&self) -> bool {
        self.base_input.calculate_pressure_derivative()
    }

    fn calculate_temperature_derivative(&self) -> bool {
        self.base_input.calculate_temperature_derivative()
    }

    fn calculate_specific_humidity_derivative(&self) -> bool {
        self.base_input.calculate_specific_humidity_derivative()
    }

    fn altitude_m(&self) -> ndarray::ArrayView1<'_, f64> {
        self.base_input.altitude_m()
    }

    fn pressure_pa(&self) -> Option<ndarray::ArrayView1<'_, f64>> {
        self.base_input.pressure_pa()
    }

    fn temperature_k(&self) -> Option<ndarray::ArrayView1<'_, f64>> {
        self.base_input.temperature_k()
    }

    fn wavelengths_nm(&self) -> Option<ndarray::ArrayView1<'_, f64>> {
        self.base_input.wavelengths_nm()
    }

    fn wavenumbers_cminv(&self) -> Option<ndarray::ArrayView1<'_, f64>> {
        self.hires_wvnum.view().into()
    }

    fn air_numberdensity_dict(&self) -> std::collections::HashMap<String, ndarray::Array1<f64>> {
        self.base_input.air_numberdensity_dict()
    }

    fn dry_air_numberdensity_dict(
        &self,
    ) -> std::collections::HashMap<String, ndarray::Array1<f64>> {
        self.base_input.dry_air_numberdensity_dict()
    }

    fn wavenumbers_cminv_left(&self) -> Option<ndarray::ArrayView1<'_, f64>> {
        self.base_input.wavenumbers_cminv_left()
    }

    fn wavenumbers_cminv_right(&self) -> Option<ndarray::ArrayView1<'_, f64>> {
        self.base_input.wavenumbers_cminv_right()
    }
}
