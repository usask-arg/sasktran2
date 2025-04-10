pub mod rayleigh;

use std::collections::HashMap;

use ndarray::*;

pub struct AtmosphereStorageOutputView<'a> {
    pub total_extinction: ArrayViewMut2<'a, f64>,
    pub ssa: ArrayViewMut2<'a, f64>,
    pub legendre: ArrayViewMut3<'a, f64>,
}

pub struct AtmosphereStorageOutputImmutView<'a> {
    pub total_extinction: ArrayView2<'a, f64>,
    pub ssa: ArrayView2<'a, f64>,
    pub legendre: ArrayView3<'a, f64>,
}

pub struct DerivMappingView<'py> {
    pub d_extinction: ArrayViewMut2<'py, f64>,
    pub d_ssa: ArrayViewMut2<'py, f64>,
    pub d_legendre: Option<ArrayViewMut3<'py, f64>>,
    pub scat_factor: Option<ArrayViewMut2<'py, f64>>,
}

pub trait StorageInputs {
    fn num_stokes(&self) -> usize;
    fn pressure_pa(&self) -> Option<ArrayView1<f64>>;
    fn temperature_k(&self) -> Option<ArrayView1<f64>>;
    fn wavelengths_nm(&self) -> Option<ArrayView1<f64>>;
    fn wavenumbers_cminv(&self) -> Option<ArrayView1<f64>>;
    fn air_numberdensity_dict(&self) -> HashMap<String, Array1<f64>>;
}

pub trait StorageOutputs {
    fn mut_view(&mut self) -> AtmosphereStorageOutputView;
    fn view(&self) -> AtmosphereStorageOutputImmutView;
}

pub trait DerivMapping<'py> {
    fn with_scatterer(self) -> Self;
    fn mut_view(&mut self) -> DerivMappingView;
    fn set_interpolator(&mut self, interpolator: &Array2<f64>);
    fn set_interp_dim(&mut self, interp_dim: &str);
    fn set_assign_name(&mut self, assign_name: &str);
}

pub trait DerivMappingGenerator<'a> {
    fn get_derivative_mapping(&self, name: &str) -> impl DerivMapping<'a>;
}

pub trait Constituent {
    fn add_to_atmosphere(
        &self,
        inputs: &impl StorageInputs,
        outputs: &mut impl StorageOutputs,
    );

    fn register_derivatives<'a>(
        &self,
        inputs: &impl StorageInputs,
        outputs: &impl StorageOutputs,
        constituent_name: &str,
        derivative_generator: &impl DerivMappingGenerator<'a>,
    );
}
