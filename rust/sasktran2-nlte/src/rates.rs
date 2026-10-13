//! Temperature-dependent rate coefficients.

/// `k(T) = a (T / t0)^n exp(-ea_over_r / T)` in SI units
/// (s^-1 m^(3(order-1))).
///
/// A constant rate is `n = 0, ea_over_r = 0`. The `(T/t0)^n` and exponential
/// factors are skipped when they are identically one.
#[derive(Clone, Debug, PartialEq)]
pub struct RateLaw {
    pub a: f64,
    pub n: f64,
    pub t0: f64,
    pub ea_over_r: f64,
}

impl RateLaw {
    pub fn constant(value: f64) -> Self {
        Self {
            a: value,
            n: 0.0,
            t0: 300.0,
            ea_over_r: 0.0,
        }
    }

    pub fn evaluate(&self, temperature_k: f64) -> f64 {
        let mut k = self.a;
        if self.n != 0.0 {
            k *= (temperature_k / self.t0).powf(self.n);
        }
        if self.ea_over_r != 0.0 {
            k *= (-self.ea_over_r / temperature_k).exp();
        }
        k
    }
}

/// Reaction order and SI conversion factor for a rate-coefficient unit string.
pub(crate) fn order_and_si_factor(units: &str) -> Option<(usize, f64)> {
    let normalised: Vec<&str> = units.split_whitespace().collect();
    match normalised.as_slice() {
        ["s-1"] => Some((1, 1.0)),
        ["m3", "s-1"] => Some((2, 1.0)),
        ["cm3", "s-1"] => Some((2, 1.0e-6)),
        ["m6", "s-1"] => Some((3, 1.0)),
        ["cm6", "s-1"] => Some((3, 1.0e-12)),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constant_rate_ignores_temperature() {
        let k = RateLaw::constant(4.0e-12);
        assert_eq!(k.evaluate(150.0), 4.0e-12);
        assert_eq!(k.evaluate(300.0), 4.0e-12);
    }

    #[test]
    fn arrhenius_with_power_law() {
        let k = RateLaw {
            a: 2.0,
            n: 0.5,
            t0: 300.0,
            ea_over_r: -67.0,
        };
        let expected = 2.0 * (200.0_f64 / 300.0).sqrt() * (67.0_f64 / 200.0).exp();
        assert!((k.evaluate(200.0) - expected).abs() < 1e-15 * expected);
    }

    #[test]
    fn unit_strings() {
        assert_eq!(order_and_si_factor("s-1"), Some((1, 1.0)));
        assert_eq!(order_and_si_factor("cm3  s-1"), Some((2, 1.0e-6)));
        assert_eq!(order_and_si_factor("cm6 s-1"), Some((3, 1.0e-12)));
        assert_eq!(order_and_si_factor("cm3/s"), None);
    }
}
