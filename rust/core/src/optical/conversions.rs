/// Converts wavelength specified in Air at STP to wavelengths at vacuum
///
/// Parameters
/// ----------
/// wavelength_nm : float
///     Wavelength in air at STP
///
/// Returns
/// -------
/// float
///    Vacuum wavelengths
///
pub fn air_wavelength_to_vacuum_wavelength(wavelength_nm: f64) -> f64 {
    let s = 1e4 / (wavelength_nm * 10.0); // convert nm to microns

    let refac_index = 1.0
        + 0.00008336624212083
        + 0.02408926869968 / (130.1065924522 - s.powf(2.0))
        + 0.0001599740894897 / (38.92568793293 - s.powf(2.0));

    wavelength_nm * refac_index
}
