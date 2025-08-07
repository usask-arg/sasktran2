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
}
