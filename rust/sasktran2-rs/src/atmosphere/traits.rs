use crate::prelude::*;

/// Mutable view of the atmosphere storage output (ssa, ext, legendre)
pub struct AtmosphereStorageOutputView<'a> {
    pub total_extinction: ArrayViewMut2<'a, f64>,
    pub ssa: ArrayViewMut2<'a, f64>,
    pub legendre: ArrayViewMut3<'a, f64>,
    pub emission_source: ArrayViewMut2<'a, f64>,
}

/// Immutable view of the atmosphere storage output (ssa, ext, legendre)
pub struct AtmosphereStorageOutputImmutView<'a> {
    pub total_extinction: ArrayView2<'a, f64>,
    pub ssa: ArrayView2<'a, f64>,
    pub legendre: ArrayView3<'a, f64>,
    pub emission_source: ArrayView2<'a, f64>,
}

/// Mutable view of the derivative mapping (d_extinction, d_ssa, d_legendre)
pub struct DerivMappingView<'py> {
    pub d_extinction: ArrayViewMut2<'py, f64>,
    pub d_ssa: ArrayViewMut2<'py, f64>,
    pub d_legendre: Option<ArrayViewMut3<'py, f64>>,
    pub scat_factor: Option<ArrayViewMut2<'py, f64>>,
    pub d_emission: ArrayViewMut2<'py, f64>,
}

/// All of the atmosphere inputs from the storage that the various constituents
/// may need to know about
pub trait StorageInputs {
    // Configuration options
    /// Number of stokes parameters
    fn num_stokes(&self) -> usize;

    /// Number of single scatter moments, note this includes the polarization factor
    /// i.e. NSTOKES=3 would have 4 * num_legendre_order moments
    fn num_singlescatter_moments(&self) -> usize;

    /// True if the user requests calculation of pressure derivatives
    fn calculate_pressure_derivative(&self) -> bool;

    /// True if the user requests calculation of temperature derivatives
    fn calculate_temperature_derivative(&self) -> bool;

    /// True if the user requests calculation of altitude derivatives
    fn calculate_specific_humidity_derivative(&self) -> bool;

    // Geometry properties
    /// Altitudes in meters of the grid
    fn altitude_m(&self) -> ArrayView1<'_, f64>;

    // Atmospheric properties
    /// Presusre in pa at altitude_m
    fn pressure_pa(&self) -> Option<ArrayView1<'_, f64>>;

    /// Temperaure in kelvin at altitude_m
    fn temperature_k(&self) -> Option<ArrayView1<'_, f64>>;

    /// Wavelengths in nm
    fn wavelengths_nm(&self) -> Option<ArrayView1<'_, f64>>;

    /// Wavenumbers in cm^-1
    fn wavenumbers_cminv(&self) -> Option<ArrayView1<'_, f64>>;

    /// hashmap of air number density factors ("N", "dN_dT", "dN_dP", "dN_dq")
    /// This is number density of air (including water vapor)
    fn air_numberdensity_dict(&self) -> HashMap<String, Array1<f64>>;

    /// hashmap of air number density factors ("N", "dN_dT", "dN_dP")
    /// This is number density of dry air
    fn dry_air_numberdensity_dict(&self) -> HashMap<String, Array1<f64>>;

    /// Allows to get parameters by string name instead of function
    fn get_parameter(&self, name: &str) -> Option<ArrayView1<'_, f64>> {
        match name {
            "pressure_pa" => self.pressure_pa(),
            "temperature_k" => self.temperature_k(),
            "wavelengths_nm" => self.wavelengths_nm(),
            "wavenumbers_cminv" => self.wavenumbers_cminv(),
            "altitude_m" => Some(self.altitude_m()),
            _ => None,
        }
    }
}

/// Trait allowing an object to return back views into the atmosphere storage
pub trait StorageOutputs {
    fn mut_view(&mut self) -> AtmosphereStorageOutputView<'_>;
    fn view(&self) -> AtmosphereStorageOutputImmutView<'_>;
}

/// Trait indicating a single derivative mapping
pub trait DerivMapping<'py> {
    fn with_scatterer(self) -> Self;
    fn mut_view(&mut self) -> DerivMappingView<'_>;
    fn set_interpolator(&mut self, interpolator: &Array2<f64>);
    fn set_interp_dim(&mut self, interp_dim: &str);
    fn set_assign_name(&mut self, assign_name: &str);
}

/// Trait for an object that is able to generate derivative mappings
pub trait DerivMappingGenerator<'a> {
    fn get_derivative_mapping(&self, name: &str) -> impl DerivMapping<'a>;
}

/// Trait for an object that is able to provide access to the atmosphere storage
pub trait AtmosphereStorageAccess {
    fn split_inputs_outputs(&mut self) -> (&impl StorageInputs, &mut impl StorageOutputs);
    fn split_inputs_outputs_deriv<'a>(
        &'a self,
    ) -> (
        &'a impl StorageInputs,
        &'a impl StorageOutputs,
        &'a impl DerivMappingGenerator<'a>,
    );
}
