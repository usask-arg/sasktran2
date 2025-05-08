use std::f64::consts::PI;

/// Bates parameterization for Rayleigh scattering
/// for the O2 refractive index
fn o2_refrac_bates(wavelength_um: f64) -> f64 {
    let coeff_ranges = [
        ([23796.7, 168988.4], (0.0, 0.221)),
        ([22120.4, 203187.6], (0.221, 0.288)),
        ([20564.8, 248089.9], (0.288, 0.546)),
        ([21351.1, 218567.0], (0.546, f64::INFINITY)),
    ];

    for (coeff, (low, high)) in coeff_ranges.iter() {
        if wavelength_um > *low && wavelength_um <= *high {
            return coeff[0] + coeff[1] / (40.9 - wavelength_um.powf(-2.0));
        }
    }

    0.0 // or maybe `panic!()` or `None` if you make this return Option<f64>
}

/// Bates parameterization for Rayleigh scattering
/// for the N2 refractive index
fn n2_refrac_bates(wavelengths_um: f64) -> f64 {
    let coeff_ranges = vec![
        ([6998.749, 3233582.0], (0.0, 0.254)),
        ([5989.242, 3363266.3], (0.254, 0.468)),
        ([6855.200, 3243157.0], (0.468, f64::INFINITY)),
    ];
    let wl = wavelengths_um;

    for (coeff, (low, high)) in coeff_ranges {
        if wl > low && wl <= high {
            let delta_lambda = 0.468 - wl;
            return coeff[0]
                + coeff[1] / (144.0 - wl.powf(-2.0))
                + 2.27684009 * delta_lambda.signum() * (-delta_lambda.abs() / 0.003).exp();
        }
    }

    0.0
}

/// Bates parameterization for Rayleigh scattering
/// for the Ar refractive index
fn ar_refrac_bates(wavelengths_um: f64) -> f64 {
    let wl = wavelengths_um;
    let nsq_m_1 = 5.547e-4 * (1.0 + 5.15e-3 * wl.powf(-2.0) + 4.19e-5 * wl.powf(-4.0));

    ((nsq_m_1 + 1.0).sqrt() - 1.0) * 1.0e8
}

/// Bates parameterization for Rayleigh scattering
/// for the CO2 refractive index
fn co2_refrac_bates(wavelength_um: f64) -> f64 {
    let wl = wavelength_um;

    22822.1
        + 117.8 * wl.powf(-2.0)
        + 2406030.0 / (130.0 - wl.powf(-2.0))
        + 15997.0 / (38.9 - wl.powf(-2.0))
}

/// Bates parameterization for Rayleigh scattering
/// for the O2 King factor
fn o2_king_bates(wavelength_um: f64) -> f64 {
    1.096 + 1.385e-3 * wavelength_um.powf(-2.0) + 1.448e-4 * wavelength_um.powf(-4.0)
}

/// Bates parameterization for Rayleigh scattering
/// for the N2 King factor
fn n2_king_bates(wavelength_um: f64) -> f64 {
    1.034 + 3.17e-4 * wavelength_um.powf(-2.0)
}

/// Bates parameterization for Rayleigh scattering
/// for the Ar King factor
fn ar_king_bates(_wavelength_um: f64) -> f64 {
    1.0
}

/// Bates parameterization for Rayleigh scattering
/// for the CO2 King factor
fn co2_king_bates(_wavelength_um: f64) -> f64 {
    1.15
}

/// Bates parameterization for Rayleigh scattering
/// Calculates (cross section, King factor) in (m^2/molecule, unitless)
/// for the given wavelength in micrometers and the
/// percentages of N2, O2, Ar, and CO2
/// in the atmosphere. The percentages should be between 0 and 100.
pub fn rayleigh_cross_section_bates(
    wavelength_um: f64,
    n2_percentage: f64,
    o2_percentage: f64,
    ar_percentage: f64,
    co2_percentage: f64,
) -> (f64, f64) {
    let lorenz_factors = o2_percentage / 100.0
        * o2_refrac_bates(wavelength_um).powf(2.0)
        * o2_king_bates(wavelength_um)
        + n2_percentage / 100.0
            * n2_refrac_bates(wavelength_um).powf(2.0)
            * n2_king_bates(wavelength_um)
        + ar_percentage / 100.0
            * ar_refrac_bates(wavelength_um).powf(2.0)
            * ar_king_bates(wavelength_um)
        + co2_percentage / 100.0
            * co2_refrac_bates(wavelength_um).powf(2.0)
            * co2_king_bates(wavelength_um);

    let eff_king = o2_percentage / 100.0 * o2_king_bates(wavelength_um)
        + n2_percentage / 100.0 * n2_king_bates(wavelength_um)
        + ar_percentage / 100.0 * ar_king_bates(wavelength_um)
        + co2_percentage / 100.0 * co2_king_bates(wavelength_um);

    let num_dens: f64 = 2.686_780_111_798_444e25;

    (
        32.0 * PI.powf(3.0) / (3.0 * num_dens.powf(2.0) * wavelength_um.powf(4.0))
            * lorenz_factors
            * 1e8,
        eff_king,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_o2_refrac_bates() {
        let result = o2_refrac_bates(0.3);
        assert!(result > 0.0, "Expected a positive refractive index for O2");
    }

    #[test]
    fn test_n2_refrac_bates() {
        let result = n2_refrac_bates(0.4);
        assert!(result > 0.0, "Expected a positive refractive index for N2");
    }

    #[test]
    fn test_ar_refrac_bates() {
        let result = ar_refrac_bates(0.5);
        assert!(result > 0.0, "Expected a positive refractive index for Ar");
    }

    #[test]
    fn test_co2_refrac_bates() {
        let result = co2_refrac_bates(0.6);
        assert!(result > 0.0, "Expected a positive refractive index for CO2");
    }

    #[test]
    fn test_o2_king_bates() {
        let result = o2_king_bates(0.3);
        assert!(result > 1.0, "Expected a King factor greater than 1 for O2");
    }

    #[test]
    fn test_n2_king_bates() {
        let result = n2_king_bates(0.4);
        assert!(result > 1.0, "Expected a King factor greater than 1 for N2");
    }

    #[test]
    fn test_ar_king_bates() {
        let result = ar_king_bates(0.5);
        assert_eq!(result, 1.0, "Expected a King factor of 1 for Ar");
    }

    #[test]
    fn test_co2_king_bates() {
        let result = co2_king_bates(0.6);
        assert_eq!(result, 1.15, "Expected a King factor of 1.15 for CO2");
    }

    #[test]
    fn test_rayleigh_cross_section_bates() {
        let (cross_section, king_factor) =
            rayleigh_cross_section_bates(0.5, 78.08, 20.95, 0.93, 0.04);
        assert!(cross_section > 0.0, "Expected a positive cross section");
        assert!(king_factor > 1.0, "Expected a King factor greater than 1");
    }
}
