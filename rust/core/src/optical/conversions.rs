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

/// Converts wavelength specified in vacuum to wavelengths in Air at STP
///
/// Parameters
/// ----------
/// wavelength_nm : float
///     Wavelength in vacuum
///
/// Returns
/// -------
/// float
///    Wavelength in air at STP
///
pub fn vacuum_wavelength_to_air_wavelength(wavelength_nm: f64) -> f64 {
    let mut air_wavelength = wavelength_nm;
    for _ in 0..10 {
        let s = 1e4 / (air_wavelength * 10.0); // convert nm to microns

        let refac_index = 1.0
            + 0.00008336624212083
            + 0.02408926869968 / (130.1065924522 - s.powf(2.0))
            + 0.0001599740894897 / (38.92568793293 - s.powf(2.0));

        air_wavelength = wavelength_nm / refac_index;
    }
    air_wavelength
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_air_wavelength_to_vacuum_wavelength() {
        let air_wavelength = 500.0; // nm
        let vacuum_wavelength = air_wavelength_to_vacuum_wavelength(air_wavelength);
        assert!((vacuum_wavelength - 500.1).abs() < 0.1); // Example tolerance
    }

    #[test]
    fn test_vacuum_wavelength_to_air_wavelength() {
        let vacuum_wavelength = 500.0; // nm
        let air_wavelength = vacuum_wavelength_to_air_wavelength(vacuum_wavelength);
        assert!((air_wavelength - 499.9).abs() < 0.1); // Example tolerance
    }

    #[test]
    fn test_round_trip_conversion() {
        let air_wavelength = 600.0; // nm
        let vacuum_wavelength = air_wavelength_to_vacuum_wavelength(air_wavelength);
        let converted_air_wavelength = vacuum_wavelength_to_air_wavelength(vacuum_wavelength);
        assert!((air_wavelength - converted_air_wavelength).abs() < 1e-6);
    }

    #[test]
    fn test_edge_case_low_wavelength() {
        let air_wavelength = 200.0; // nm
        let vacuum_wavelength = air_wavelength_to_vacuum_wavelength(air_wavelength);
        assert!(vacuum_wavelength > air_wavelength);
    }

    #[test]
    fn test_edge_case_high_wavelength() {
        let air_wavelength = 2000.0; // nm
        let vacuum_wavelength = air_wavelength_to_vacuum_wavelength(air_wavelength);
        assert!(vacuum_wavelength > air_wavelength);
    }
}
