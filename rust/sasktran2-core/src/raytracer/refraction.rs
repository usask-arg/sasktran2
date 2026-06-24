use std::error::Error;
use std::fmt;

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
        self.integrate_path_with_tangent_index(
            tangent_radius,
            self.refractive_index_at_radius(tangent_radius),
            r1,
            r2,
        )
    }

    pub(crate) fn integrate_path_with_tangent_index(
        &self,
        tangent_radius: f64,
        tangent_refractive_index: f64,
        r1: f64,
        r2: f64,
    ) -> PathIntegral {
        integrate_path(self, tangent_radius, tangent_refractive_index, r1, r2)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PathIntegral {
    pub path_length: f64,
    pub deflection_angle: f64,
}

#[allow(clippy::excessive_precision)]
const GQ64_NODES: [f64; 32] = [
    0.99930504173577213946,
    0.99634011677195527935,
    0.99101337147674432074,
    0.98333625388462595693,
    0.97332682778991096374,
    0.96100879965205371892,
    0.94641137485840281606,
    0.92956917213193957582,
    0.91052213707850280576,
    0.88931544599511410585,
    0.86599939815409281976,
    0.84062929625258036275,
    0.81326531512279755974,
    0.78397235894334140761,
    0.75281990726053189661,
    0.71988185017161082685,
    0.68523631305423324256,
    0.64896547125465733986,
    0.61115535517239325025,
    0.57189564620263403428,
    0.53127946401989454566,
    0.48940314570705295748,
    0.44636601725346408798,
    0.40227015796399160370,
    0.35722015833766811595,
    0.31132287199021095616,
    0.26468716220876741637,
    0.21742364374000708415,
    0.16964442042399281804,
    0.12146281929612055447,
    0.07299312178779903945,
    0.024350292663424432509,
];

#[allow(clippy::excessive_precision)]
const GQ64_WEIGHTS: [f64; 32] = [
    0.0017832807216964329473,
    0.0041470332605624676353,
    0.0065044579689783628561,
    0.0088467598263639477231,
    0.011168139460131128819,
    0.013463047896718642598,
    0.015726030476024719322,
    0.017951715775697343085,
    0.020134823153530209372,
    0.022270173808383254159,
    0.024352702568710873338,
    0.026377469715054658672,
    0.028339672614259483228,
    0.030234657072402478868,
    0.032057928354851553585,
    0.033805161837141609392,
    0.035472213256882383811,
    0.03705512854024004604,
    0.038550153178615629129,
    0.039953741132720341387,
    0.04126256324262352861,
    0.042473515123653589007,
    0.043583724529323453377,
    0.04459055816375656306,
    0.04549162792741814448,
    0.046284796581314417296,
    0.046968182816210017325,
    0.047540165714830308662,
    0.047999388596458307728,
    0.04834476223480295717,
    0.048575467441503426935,
    0.048690957009139720383,
];

pub fn integrate_path(
    profile: &RefractiveProfile,
    tangent_radius: f64,
    nt: f64,
    mut r1: f64,
    mut r2: f64,
) -> PathIntegral {
    const MIN_CELL_LENGTH: f64 = 0.1;
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

    let mut extra_path = 0.0;
    let mut deflection_angle = 0.0;
    for (&mu, &weight) in GQ64_NODES.iter().zip(GQ64_WEIGHTS.iter()) {
        let a1 = 0.5 * mu + 0.5;
        let a2 = -0.5 * mu + 0.5;
        let a3 = 0.5 * mu - 0.5;
        let a4 = -0.5 * mu - 0.5;
        let w = 0.5 * weight;

        let x1 = half_width * a1 + center;
        let x2 = half_width * a2 + center;
        let x3 = half_width * a3 + center;
        let x4 = half_width * a4 + center;

        let n1 = profile.refractive_index_at_radius(x1 * x1 + tangent_radius);
        let n2 = profile.refractive_index_at_radius(x2 * x2 + tangent_radius);
        let n3 = profile.refractive_index_at_radius(x3 * x3 + tangent_radius);
        let n4 = profile.refractive_index_at_radius(x4 * x4 + tangent_radius);

        extra_path += w * path_integrand(x1, n1);
        extra_path += w * path_integrand(x2, n2);
        extra_path += w * path_integrand(x3, n3);
        extra_path += w * path_integrand(x4, n4);

        deflection_angle += w * angle_integrand(x1, n1);
        deflection_angle += w * angle_integrand(x2, n2);
        deflection_angle += w * angle_integrand(x3, n3);
        deflection_angle += w * angle_integrand(x4, n4);
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

    #[test]
    fn nontrivial_profile_path_integral_is_finite() {
        let profile =
            RefractiveProfile::new(10.0, vec![0.0, 10.0, 20.0], vec![1.1, 1.05, 1.0]).unwrap();
        let result = profile.integrate_path(12.0, 15.0, 25.0);

        assert!(result.path_length.is_finite());
        assert!(result.deflection_angle.is_finite());
        assert!(result.path_length > 0.0);
        assert!(result.deflection_angle > 0.0);
    }

    #[test]
    fn near_tangent_integral_applies_radius_dither() {
        let profile =
            RefractiveProfile::new(10.0, vec![0.0, 10.0, 20.0], vec![1.1, 1.05, 1.0]).unwrap();
        let result = profile.integrate_path(15.0, 15.0, 15.05);

        assert!(result.path_length.is_finite());
        assert!(result.deflection_angle.is_finite());
        assert!(result.path_length >= 0.0);
    }
}
