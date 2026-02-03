use crate::optical::storage::OpticalQuantities;
use anyhow::*;
use gauss_quad::GaussLegendre;
use ndarray::{Array1, ArrayView1};
use rebasis::grid::Grid;
use crate::util::{argsort_f64, is_sorted};

use super::StorageInputs;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LineshapeCoordinate {
    WavenumberCminv,
    WavelengthNm,
}

/// The internal spectral grid object used for atmosphere storage.
/// The user can specify either wavenumber or wavelength, but internally wavenumbers is preferred
/// Each wavenumber point can either represent a delta function (monochromatic) or an integrated lineshape
pub struct SpectralGrid {
    central_wavenumbers_cminv: Array1<f64>,
    central_wavelengths_nm: Array1<f64>,
    left_wavenumbers_cminv: Array1<f64>,
    right_wavenumbers_cminv: Array1<f64>,
    left_wavelengths_nm: Array1<f64>,
    right_wavelengths_nm: Array1<f64>,
    basis_grid: Grid,
    wavenumber_sort_idx: Vec<usize>,
    sorted: bool,
    coordinate: LineshapeCoordinate,
}

impl SpectralGrid {
    pub fn from_grid(grid: Grid, natural_coordinate: LineshapeCoordinate) -> Result<Self> {
        match natural_coordinate {
            LineshapeCoordinate::WavenumberCminv => {
                let center_wave = Array1::from_vec(grid.vec_basis().iter().map(|b| {
                    b.center()
                }).collect());

                let left_wave = Array1::from_vec(grid.vec_basis().iter().map(|b| {
                    b.lower_limit()
                }).collect());

                let right_wave = Array1::from_vec(grid.vec_basis().iter().map(|b| {
                    b.upper_limit()
                }).collect());

                let left_wavelength = right_wave.mapv(|x| 1.0e7 / x);
                let right_wavelength = left_wave.mapv(|x| 1.0e7 / x);

                let center_wavelength = center_wave.mapv(|x| 1.0e7 / x);

                let sort_idx = argsort_f64(center_wave.as_slice().unwrap());
                let sorted = is_sorted(center_wave.as_slice().unwrap());

                Ok(SpectralGrid {
                    central_wavenumbers_cminv: center_wave,
                    central_wavelengths_nm: center_wavelength,
                    left_wavenumbers_cminv: left_wave,
                    right_wavenumbers_cminv: right_wave,
                    left_wavelengths_nm: left_wavelength,
                    right_wavelengths_nm: right_wavelength,
                    basis_grid: grid,
                    wavenumber_sort_idx: sort_idx,
                    sorted,
                    coordinate: LineshapeCoordinate::WavenumberCminv,
                })
            }
            LineshapeCoordinate::WavelengthNm => {
                let center_wave = Array1::from_vec(grid.vec_basis().iter().map(|b| {
                    b.center()
                }).collect());

                let left_wave = Array1::from_vec(grid.vec_basis().iter().map(|b| {
                    b.lower_limit()
                }).collect());

                let right_wave = Array1::from_vec(grid.vec_basis().iter().map(|b| {
                    b.upper_limit()
                }).collect());

                let left_wavenumber = right_wave.mapv(|x| 1.0e7 / x);
                let right_wavenumber = left_wave.mapv(|x| 1.0e7 / x);

                let center_wavenumber = center_wave.mapv(|x| 1.0e7 / x);

                let sort_idx = argsort_f64(center_wave.as_slice().unwrap());
                let sorted = is_sorted(center_wave.as_slice().unwrap());

                Ok(SpectralGrid {
                    central_wavenumbers_cminv: center_wavenumber,
                    central_wavelengths_nm: center_wave,
                    left_wavenumbers_cminv: left_wavenumber,
                    right_wavenumbers_cminv: right_wavenumber,
                    left_wavelengths_nm: left_wave,
                    right_wavelengths_nm: right_wave,
                    basis_grid: grid,
                    wavenumber_sort_idx: sort_idx,
                    sorted,
                    coordinate: LineshapeCoordinate::WavelengthNm,
                })
            }
        }


    }

    pub fn from_wavelengths_nm(wavelengths_nm: &Array1<f64>) -> Result<Self> {
        let basis_vec = wavelengths_nm
            .iter()
            .map(|&wvnum| rebasis::basis::Delta::new(wvnum ))
            .map(rebasis::basis::BasisType::Delta)
            .collect();

        let basis_grid = rebasis::grid::Grid::new(basis_vec);

        SpectralGrid::from_grid(basis_grid, LineshapeCoordinate::WavelengthNm)
    }

    pub fn from_wavenumbers_cminv(wavenumbers_cminv: &Array1<f64>) -> Result<Self> {
        let basis_vec = wavenumbers_cminv
            .iter()
            .map(|&wvnum| rebasis::basis::Delta::new(wvnum ))
            .map(rebasis::basis::BasisType::Delta)
            .collect();

        let basis_grid = rebasis::grid::Grid::new(basis_vec);

        SpectralGrid::from_grid(basis_grid, LineshapeCoordinate::WavenumberCminv)
    }

    pub fn central_wavenumber_cminv(&self) -> ArrayView1<'_, f64> {
        self.central_wavenumbers_cminv.view()
    }

    pub fn central_wavelengths_nm(&self) -> ArrayView1<'_, f64> {
        self.central_wavelengths_nm.view()
    }

    pub fn left_wavenumbers_cminv(&self) -> ArrayView1<'_, f64> {
        self.left_wavenumbers_cminv.view()
    }

    pub fn right_wavenumbers_cminv(&self) -> ArrayView1<'_, f64> {
        self.right_wavenumbers_cminv.view()
    }

    pub fn left_wavelengths_nm(&self) -> ArrayView1<'_, f64> {
        self.left_wavelengths_nm.view()
    }

    pub fn right_wavelengths_nm(&self) -> ArrayView1<'_, f64> {
        self.right_wavelengths_nm.view()
    }

    pub fn basis_grid(&self) -> &Grid {
        &self.basis_grid
    }

    pub fn sort_idx(&self) -> &Vec<usize> {
        &self.wavenumber_sort_idx
    }

    pub fn coordinate(&self) -> LineshapeCoordinate {
        self.coordinate
    }
}

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
    spectral_grid: Option<SpectralGrid>,
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
        self.spectral_grid = Some(SpectralGrid::from_wavelengths_nm(&wavelengths_nm).unwrap());

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

    fn air_numberdensity_dict(&self) -> std::collections::HashMap<String, ndarray::Array1<f64>> {
        todo!()
    }

    fn dry_air_numberdensity_dict(
        &self,
    ) -> std::collections::HashMap<String, ndarray::Array1<f64>> {
        todo!()
    }

    fn spectral_grid(&self) -> Option<&SpectralGrid> {
        self.spectral_grid.as_ref()
    }

    fn spectral_integration_mode(&self) -> crate::bindings::config::SpectralGridMode {
        crate::bindings::config::SpectralGridMode::Monochromatic
    }

}
