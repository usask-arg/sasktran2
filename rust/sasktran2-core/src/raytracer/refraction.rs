use std::error::Error;
use std::fmt;

use gauss_quad::legendre::GaussLegendre;

#[derive(Debug, Clone, PartialEq)]
pub enum RefractiveProfileError {
    MismatchedLength {
        altitudes: usize,
        refractive_index: usize,
    },
    TooFewPoints,
    NonIncreasingAltitude {
        index: usize,
    },
    NonPositiveRefractiveIndex {
        index: usize,
    },
}

impl fmt::Display for RefractiveProfileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MismatchedLength {
                altitudes,
                refractive_index,
            } => write!(
                f,
                "altitude length {altitudes} does not match refractive-index length {refractive_index}"
            ),
            Self::TooFewPoints => write!(f, "a refractive profile needs at least two points"),
            Self::NonIncreasingAltitude { index } => {
                write!(
                    f,
                    "altitudes must be strictly increasing near index {index}"
                )
            }
            Self::NonPositiveRefractiveIndex { index } => {
                write!(f, "refractive index must be positive at index {index}")
            }
        }
    }
}

impl Error for RefractiveProfileError {}

#[derive(Debug, Clone, PartialEq)]
pub struct RefractiveProfile {
    earth_radius: f64,
    altitudes: Vec<f64>,
    refractive_index: Vec<f64>,
}

impl RefractiveProfile {
    pub fn new(
        earth_radius: f64,
        altitudes: Vec<f64>,
        refractive_index: Vec<f64>,
    ) -> Result<Self, RefractiveProfileError> {
        if altitudes.len() != refractive_index.len() {
            return Err(RefractiveProfileError::MismatchedLength {
                altitudes: altitudes.len(),
                refractive_index: refractive_index.len(),
            });
        }
        if altitudes.len() < 2 {
            return Err(RefractiveProfileError::TooFewPoints);
        }
        for i in 1..altitudes.len() {
            if altitudes[i] <= altitudes[i - 1] {
                return Err(RefractiveProfileError::NonIncreasingAltitude { index: i });
            }
        }
        for (index, value) in refractive_index.iter().enumerate() {
            if *value <= 0.0 {
                return Err(RefractiveProfileError::NonPositiveRefractiveIndex { index });
            }
        }
        Ok(Self {
            earth_radius,
            altitudes,
            refractive_index,
        })
    }

    #[inline(always)]
    pub fn earth_radius(&self) -> f64 {
        self.earth_radius
    }

    #[inline(always)]
    pub fn altitudes(&self) -> &[f64] {
        &self.altitudes
    }

    #[inline(always)]
    pub fn refractive_index(&self) -> &[f64] {
        &self.refractive_index
    }

    #[inline]
    pub fn refractive_index_at_radius(&self, radius: f64) -> f64 {
        self.refractive_index_at_altitude(radius - self.earth_radius)
    }

    pub fn refractive_index_at_altitude(&self, altitude: f64) -> f64 {
        if altitude <= self.altitudes[0] {
            return self.refractive_index[0];
        }
        let last = self.altitudes.len() - 1;
        if altitude >= self.altitudes[last] {
            return self.refractive_index[last];
        }

        let upper = self
            .altitudes
            .partition_point(|boundary| *boundary <= altitude);
        let lower = upper - 1;
        let width = self.altitudes[upper] - self.altitudes[lower];
        let upper_weight = (altitude - self.altitudes[lower]) / width;
        let log_n = (1.0 - upper_weight) * self.refractive_index[lower].ln()
            + upper_weight * self.refractive_index[upper].ln();
        log_n.exp()
    }

    pub fn tangent_radius(&self, straight_line_tangent_radius: f64) -> f64 {
        const MAX_ITERATIONS: usize = 500;
        const TOLERANCE: f64 = 1e-6;

        let mut current = straight_line_tangent_radius;
        let mut next = current;

        for _ in 0..MAX_ITERATIONS {
            let n = self.refractive_index_at_radius(current);
            next = straight_line_tangent_radius / n;

            if (next - current).abs() < TOLERANCE {
                return current;
            }
            current = next;
        }

        if (next - current).abs() < 1.0 {
            (next + current) / 2.0
        } else {
            current
        }
    }

    pub fn integrate_path(&self, tangent_radius: f64, r1: f64, r2: f64) -> PathIntegral {
        integrate_path(self, tangent_radius, r1, r2)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PathIntegral {
    pub path_length: f64,
    pub deflection_angle: f64,
}

pub fn integrate_path(
    profile: &RefractiveProfile,
    tangent_radius: f64,
    mut r1: f64,
    mut r2: f64,
) -> PathIntegral {
    const MIN_CELL_LENGTH: f64 = 0.1;
    const NUM_INTEGRATION_POINTS: usize = 64;
    const RADIUS_DITHER: f64 = 1e-6;

    if r2 < r1 {
        std::mem::swap(&mut r1, &mut r2);
    }

    let min_radius = tangent_radius + RADIUS_DITHER;
    if r1 <= min_radius {
        r1 = min_radius;
    }
    if r2 <= min_radius {
        r2 = min_radius;
    }

    if (r1 - r2).abs() < MIN_CELL_LENGTH {
        return PathIntegral {
            path_length: (r2 * r2 - tangent_radius * tangent_radius).sqrt()
                - (r1 * r1 - tangent_radius * tangent_radius).sqrt(),
            deflection_angle: (tangent_radius / r2).acos() - (tangent_radius / r1).acos(),
        };
    }

    let nt = profile.refractive_index_at_radius(tangent_radius);

    let path_integrand = |x: f64, n: f64| {
        let sqf = (1.0 + (n - nt) / n * tangent_radius / (x * x)).sqrt()
            * (x * x + (n + nt) / n * tangent_radius).sqrt();
        let f = (x * x + 2.0 * tangent_radius).sqrt() + sqf;
        let g = (x * x + 2.0 * tangent_radius).sqrt() * sqf;

        2.0 * tangent_radius * tangent_radius * (nt + n) / n * (nt - n) / n
            * (x * x + tangent_radius)
            / (x * x * f * g)
    };

    let angle_integrand = |x: f64, n: f64| {
        let r = x * x + tangent_radius;
        let t1 = ((n - nt) / n * r / (x * x) + nt / n).sqrt();
        let t2 = (r + tangent_radius * nt / n).sqrt();

        2.0 * nt * tangent_radius / r * (1.0 / (n * t1 * t2))
    };

    let x_low = (r1 - tangent_radius).sqrt();
    let x_high = (r2 - tangent_radius).sqrt();
    let half_width = (x_high - x_low) / 2.0;
    let center = (x_high + x_low) / 2.0;
    let quad = GaussLegendre::new(NUM_INTEGRATION_POINTS)
        .expect("64-point Gauss-Legendre quadrature should be constructible");

    let mut extra_path = 0.0;
    let mut deflection_angle = 0.0;
    for (node, weight) in quad.iter() {
        let x = half_width * *node + center;
        let n = profile.refractive_index_at_radius(x * x + tangent_radius);

        extra_path += *weight * path_integrand(x, n);
        deflection_angle += *weight * angle_integrand(x, n);
    }

    extra_path *= half_width;
    deflection_angle *= half_width;

    let straight_path = (r2 * r2 - tangent_radius * tangent_radius).sqrt()
        - (r1 * r1 - tangent_radius * tangent_radius).sqrt();

    let mut result = PathIntegral {
        path_length: straight_path + extra_path,
        deflection_angle,
    };

    if !result.path_length.is_finite() || !result.deflection_angle.is_finite() {
        result.path_length = 0.0;
        result.deflection_angle = 0.0;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unity_profile() -> RefractiveProfile {
        RefractiveProfile::new(10.0, vec![0.0, 10.0, 20.0], vec![1.0, 1.0, 1.0]).unwrap()
    }

    #[test]
    fn refractive_index_interpolates_log_linearly() {
        let profile = RefractiveProfile::new(0.0, vec![0.0, 10.0], vec![1.0, 4.0]).unwrap();

        assert!((profile.refractive_index_at_altitude(5.0) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn tangent_radius_is_unchanged_for_unity_index() {
        let profile = unity_profile();
        assert!((profile.tangent_radius(12.0) - 12.0).abs() < 1e-12);
    }

    #[test]
    fn tangent_radius_decreases_for_index_above_one() {
        let profile =
            RefractiveProfile::new(10.0, vec![0.0, 10.0, 20.0], vec![1.1, 1.05, 1.0]).unwrap();

        assert!(profile.tangent_radius(20.0) < 20.0);
    }

    #[test]
    fn unity_index_path_integral_matches_straight_path() {
        let profile = unity_profile();
        let result = profile.integrate_path(10.0, 15.0, 20.0);
        let expected =
            (20.0_f64 * 20.0 - 10.0 * 10.0).sqrt() - (15.0_f64 * 15.0 - 10.0 * 10.0).sqrt();

        assert!((result.path_length - expected).abs() < 1e-8);
    }
}
