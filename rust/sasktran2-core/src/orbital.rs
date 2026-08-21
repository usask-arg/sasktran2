//! Pure orbital-plane grouping, projection, and refractivity utilities.

use std::collections::BTreeMap;

use nalgebra::{Matrix3, SymmetricEigen, Vector3};

use crate::raytracer::Vec3;

const PI: f64 = std::f64::consts::PI;
const TWO_PI: f64 = 2.0 * PI;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrackCoordinate {
    pub segment: usize,
    pub fraction: f64,
    pub angle: f64,
}

#[derive(Debug, Clone)]
pub struct FittedPlane {
    pub reference_z: Vec3,
    pub reference_x: Vec3,
    pub horizontal_angles: Vec<f64>,
    pub max_out_of_plane_angle: f64,
}

fn normalize(value: Vec3, name: &str) -> Result<Vec3, String> {
    let norm = value.norm();
    if !value.is_finite() || !norm.is_finite() || norm == 0.0 {
        return Err(format!("{name} must be finite and non-zero"));
    }
    Ok(value / norm)
}

fn angular_separation(left: Vec3, right: Vec3) -> f64 {
    left.cross(right)
        .norm()
        .atan2(left.dot(right).clamp(-1.0, 1.0))
}

pub fn cumulative_track_angles(track: &[Vec3]) -> Result<Vec<f64>, String> {
    if track.len() < 2 {
        return Err("An orbital track requires at least two positions".to_string());
    }
    let track = track
        .iter()
        .map(|value| normalize(*value, "Orbital track direction"))
        .collect::<Result<Vec<_>, _>>()?;
    let mut cumulative = vec![0.0; track.len()];
    for index in 1..track.len() {
        let segment = track[index - 1].dot(track[index]).clamp(-1.0, 1.0).acos();
        if !segment.is_finite() || segment <= 1e-12 || segment >= PI {
            return Err(format!(
                "Orbital track segment {index} must have angular length in (0, pi)"
            ));
        }
        cumulative[index] = cumulative[index - 1] + segment;
    }
    Ok(cumulative)
}

pub fn time_groups(times_ns: &[i64], duration_ns: i64) -> Result<Vec<Vec<usize>>, String> {
    if times_ns.is_empty() {
        return Err("At least one observation is required".to_string());
    }
    if duration_ns <= 0 {
        return Err("The time-group duration must be positive".to_string());
    }
    if times_ns.contains(&i64::MIN) {
        return Err("Observation times must not contain NaT".to_string());
    }
    let first_time = *times_ns.iter().min().unwrap();
    let mut grouped = BTreeMap::<i128, Vec<usize>>::new();
    for (index, time) in times_ns.iter().enumerate() {
        let bin = (*time as i128 - first_time as i128).div_euclid(duration_ns as i128);
        grouped.entry(bin).or_default().push(index);
    }
    Ok(grouped.into_values().collect())
}

pub fn locate_track_coordinate(
    target: Vec3,
    track: &[Vec3],
    cumulative: &[f64],
) -> TrackCoordinate {
    locate_track_coordinate_near(target, track, cumulative, None)
}

/// Locate a direction on an ordered spherical polyline, using an expected
/// cumulative coordinate to disambiguate closed or self-crossing tracks.
///
/// Spatial distance remains the primary criterion. The expected coordinate is
/// only used when two candidates are coincident to numerical precision.
pub fn locate_track_coordinate_near(
    target: Vec3,
    track: &[Vec3],
    cumulative: &[f64],
    expected_angle: Option<f64>,
) -> TrackCoordinate {
    let mut best = TrackCoordinate {
        segment: 0,
        fraction: 0.0,
        angle: cumulative[0],
    };
    let mut best_residual = f64::INFINITY;
    let mut best_expected_distance = f64::INFINITY;
    for segment in 0..track.len() - 1 {
        let length = cumulative[segment + 1] - cumulative[segment];
        let tangent = (track[segment + 1] - length.cos() * track[segment]) / length.sin();
        let along = target.dot(tangent).atan2(target.dot(track[segment]));
        let clamped = along.clamp(0.0, length);
        let candidate = clamped.cos() * track[segment] + clamped.sin() * tangent;
        let residual = angular_separation(target, candidate);
        let angle = cumulative[segment] + clamped;
        let expected_distance = expected_angle.map_or(0.0, |expected| (angle - expected).abs());
        if residual < best_residual - 1e-10
            || ((residual - best_residual).abs() <= 1e-10
                && expected_distance < best_expected_distance)
        {
            best_residual = residual;
            best_expected_distance = expected_distance;
            best = TrackCoordinate {
                segment,
                fraction: clamped / length,
                angle,
            };
        }
    }
    best
}

/// Locate a direction on only the portion of an ordered spherical polyline
/// contained in the supplied cumulative-coordinate interval.
pub fn locate_track_coordinate_in_range(
    target: Vec3,
    track: &[Vec3],
    cumulative: &[f64],
    minimum_angle: f64,
    maximum_angle: f64,
) -> Result<TrackCoordinate, String> {
    locate_track_coordinate_in_range_near(
        target,
        track,
        cumulative,
        minimum_angle,
        maximum_angle,
        None,
    )
}

/// Locate a direction in a cumulative-coordinate interval and use an expected
/// coordinate to break spatially coincident ties. This is useful when mapping
/// an ordered sequence across a closed or self-crossing track.
pub fn locate_track_coordinate_in_range_near(
    target: Vec3,
    track: &[Vec3],
    cumulative: &[f64],
    minimum_angle: f64,
    maximum_angle: f64,
    expected_angle: Option<f64>,
) -> Result<TrackCoordinate, String> {
    if !minimum_angle.is_finite()
        || !maximum_angle.is_finite()
        || minimum_angle > maximum_angle
        || minimum_angle < cumulative[0]
        || maximum_angle > *cumulative.last().unwrap()
        || expected_angle.is_some_and(|expected| {
            !expected.is_finite() || expected < minimum_angle || expected > maximum_angle
        })
    {
        return Err("The track-coordinate search interval is invalid".to_string());
    }
    let mut best = None;
    let mut best_residual = f64::INFINITY;
    let mut best_expected_distance = f64::INFINITY;
    for segment in 0..track.len() - 1 {
        let segment_start = cumulative[segment].max(minimum_angle);
        let segment_end = cumulative[segment + 1].min(maximum_angle);
        if segment_start > segment_end {
            continue;
        }
        let length = cumulative[segment + 1] - cumulative[segment];
        let tangent = (track[segment + 1] - length.cos() * track[segment]) / length.sin();
        let along = target.dot(tangent).atan2(target.dot(track[segment]));
        let minimum_offset = segment_start - cumulative[segment];
        let maximum_offset = segment_end - cumulative[segment];
        let clamped = along.clamp(minimum_offset, maximum_offset);
        let candidate = clamped.cos() * track[segment] + clamped.sin() * tangent;
        let residual = angular_separation(target, candidate);
        let angle = cumulative[segment] + clamped;
        let expected_distance = expected_angle.map_or(0.0, |expected| (angle - expected).abs());
        let is_better = if expected_angle.is_some() {
            residual < best_residual - 1e-10
                || ((residual - best_residual).abs() <= 1e-10
                    && expected_distance < best_expected_distance)
        } else {
            residual < best_residual
        };
        if is_better {
            best_residual = residual;
            best_expected_distance = expected_distance;
            best = Some(TrackCoordinate {
                segment,
                fraction: clamped / length,
                angle,
            });
        }
    }
    best.ok_or_else(|| "The track-coordinate search interval contains no segment".to_string())
}

pub fn fit_great_circle_plane(
    track: &[Vec3],
    grid_indices: &[usize],
) -> Result<FittedPlane, String> {
    if grid_indices.len() < 2 {
        return Err("A projected orbital window needs at least two points".to_string());
    }
    let vectors = track
        .iter()
        .map(|value| Vector3::new(value.x, value.y, value.z))
        .collect::<Vec<_>>();
    let mut covariance = Matrix3::zeros();
    for &index in grid_indices {
        covariance += vectors[index] * vectors[index].transpose();
    }
    let eigensystem = SymmetricEigen::new(covariance);
    let normal_index = eigensystem
        .eigenvalues
        .iter()
        .enumerate()
        .min_by(|(_, left), (_, right)| left.total_cmp(right))
        .map(|(index, _)| index)
        .ok_or_else(|| "Unable to fit orbital plane".to_string())?;
    let normal = eigensystem.eigenvectors.column(normal_index).into_owned();
    let center = grid_indices.len() / 2;
    let center_direction = vectors[grid_indices[center]];
    let reference_z = (center_direction - center_direction.dot(&normal) * normal).normalize();
    let first = vectors[*grid_indices.first().unwrap()];
    let last = vectors[*grid_indices.last().unwrap()];
    let forward = last - first;
    let reference_x =
        (forward - forward.dot(&normal) * normal - forward.dot(&reference_z) * reference_z)
            .normalize();
    if !reference_z.iter().all(|value| value.is_finite())
        || !reference_x.iter().all(|value| value.is_finite())
    {
        return Err("The projected orbital-plane basis is degenerate".to_string());
    }
    let reference_y = reference_z.cross(&reference_x).normalize();

    let mut horizontal_angles = Vec::with_capacity(grid_indices.len());
    let mut previous = 0.0;
    for (local_index, &global_index) in grid_indices.iter().enumerate() {
        let direction = vectors[global_index];
        let projected = (direction - direction.dot(&reference_y) * reference_y).normalize();
        let mut angle = projected
            .dot(&reference_x)
            .atan2(projected.dot(&reference_z));
        if local_index > 0 {
            while angle <= previous - PI {
                angle += TWO_PI;
            }
            while angle > previous + PI {
                angle -= TWO_PI;
            }
            if angle <= previous + 1e-12 {
                return Err("The projected orbital window folds or is non-monotonic".to_string());
            }
        }
        previous = angle;
        horizontal_angles.push(angle);
    }
    let center_angle = horizontal_angles[center];
    for angle in &mut horizontal_angles {
        *angle -= center_angle;
    }
    if horizontal_angles.last().unwrap() - horizontal_angles.first().unwrap() >= PI {
        return Err("The projected orbital window spans pi radians or more".to_string());
    }
    let max_out_of_plane_angle = grid_indices
        .iter()
        .map(|&index| {
            vectors[index]
                .dot(&reference_y)
                .abs()
                .clamp(0.0, 1.0)
                .asin()
        })
        .fold(0.0, f64::max);
    Ok(FittedPlane {
        reference_z: Vec3::new(reference_z.x, reference_z.y, reference_z.z),
        reference_x: Vec3::new(reference_x.x, reference_x.y, reference_x.z),
        horizontal_angles,
        max_out_of_plane_angle,
    })
}

pub fn ciddor_index(
    temperature_k: f64,
    pressure_pa: f64,
    specific_humidity: f64,
    co2_ppm: f64,
    wavelength_nm: f64,
) -> Result<f64, String> {
    if !temperature_k.is_finite()
        || temperature_k <= 0.0
        || !pressure_pa.is_finite()
        || pressure_pa < 0.0
        || !specific_humidity.is_finite()
        || !(0.0..1.0).contains(&specific_humidity)
        || !co2_ppm.is_finite()
        || co2_ppm < 0.0
        || !wavelength_nm.is_finite()
        || wavelength_nm <= 0.0
    {
        return Err("Pressure, temperature, humidity, CO2, and wavelength are invalid".to_string());
    }
    let temperature_c = temperature_k - 273.15;
    let x_h2o =
        specific_humidity / (specific_humidity + (1.0 - specific_humidity) / (18.01528 / 28.9647));
    let spectral = 1.0 / (wavelength_nm / 1e3).powi(2);
    let ras = 1e-8 * (5_792_105.0 / (238.0185 - spectral) + 167_917.0 / (57.362 - spectral));
    let rvs = 1.022e-8
        * (295.235 + 2.6422 * spectral - 0.03238 * spectral.powi(2) + 0.004028 * spectral.powi(3));
    let molar_mass_air = 0.0289635 + 1.2011e-8 * (co2_ppm - 400.0);
    let raxs = ras * (1.0 + 5.34e-7 * (co2_ppm - 450.0));
    let z = 1.0
        - (pressure_pa / temperature_k)
            * (1.58123e-6 - 2.9331e-8 * temperature_c
                + 1.1043e-10 * temperature_c.powi(2)
                + (5.707e-6 - 2.051e-8 * temperature_c) * x_h2o
                + (1.9898e-4 - 2.376e-6 * temperature_c) * x_h2o.powi(2))
        + (pressure_pa / temperature_k).powi(2) * (1.83e-11 - 0.765e-8 * x_h2o.powi(2));
    let rho_axs = 101_325.0 * molar_mass_air / (0.9995922115 * 8.314472 * 288.15);
    let rho_v = x_h2o * pressure_pa * 0.018015 / (z * 8.314472 * temperature_k);
    let rho_a = (1.0 - x_h2o) * pressure_pa * molar_mass_air / (z * 8.314472 * temperature_k);
    Ok(1.0 + (rho_a / rho_axs) * raxs + (rho_v / 0.00985938) * rvs)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn great_circle() -> Vec<Vec3> {
        (0..81)
            .map(|index| {
                let angle = -1.0 + index as f64 * 0.025;
                Vec3::new(angle.sin(), 0.0, angle.cos())
            })
            .collect()
    }

    #[test]
    fn bins_are_half_open_and_anchored_at_the_earliest_sample() {
        assert_eq!(
            time_groups(&[20, 0, 10, 9], 10).unwrap(),
            vec![vec![1, 3], vec![2], vec![0]]
        );
    }

    #[test]
    fn time_groups_reject_nat_and_handle_the_full_integer_range() {
        assert!(time_groups(&[0, i64::MIN], 10).is_err());
        assert_eq!(
            time_groups(&[i64::MIN + 1, i64::MAX], i64::MAX).unwrap(),
            vec![vec![0], vec![1]]
        );
    }

    #[test]
    fn cumulative_and_tangent_coordinates_follow_track_order() {
        let track = great_circle();
        let cumulative = cumulative_track_angles(&track).unwrap();
        let angle: f64 = 0.1375;
        let target = Vec3::new(angle.sin(), 0.0, angle.cos());
        let coordinate = locate_track_coordinate(target, &track, &cumulative);
        assert!((coordinate.angle - (angle + 1.0)).abs() < 1e-12);
        assert!((coordinate.fraction - 0.5).abs() < 1e-12);
    }

    #[test]
    fn preferred_coordinate_disambiguates_a_closed_track() {
        let track = (0..=360)
            .map(|index| {
                let angle = index as f64 * TWO_PI / 360.0;
                Vec3::new(angle.sin(), 0.0, angle.cos())
            })
            .collect::<Vec<_>>();
        let cumulative = cumulative_track_angles(&track).unwrap();
        let target = Vec3::new(0.0, 0.0, 1.0);

        let first = locate_track_coordinate(target, &track, &cumulative);
        let last = locate_track_coordinate_near(target, &track, &cumulative, Some(TWO_PI));
        let ranged = locate_track_coordinate_in_range(
            target,
            &track,
            &cumulative,
            PI,
            *cumulative.last().unwrap(),
        )
        .unwrap();

        assert!(first.angle.abs() < 1e-12);
        let total = *cumulative.last().unwrap();
        assert!((total - TWO_PI).abs() < 2e-12);
        assert!((last.angle - total).abs() < 1e-12);
        assert!((ranged.angle - total).abs() < 1e-12);
    }

    #[test]
    fn ranged_preference_supports_reverse_order_on_a_closed_track() {
        let track = (0..=72)
            .map(|index| {
                let angle = index as f64 * TWO_PI / 72.0;
                Vec3::new(angle.sin(), 0.0, angle.cos())
            })
            .collect::<Vec<_>>();
        let cumulative = cumulative_track_angles(&track).unwrap();
        let total = *cumulative.last().unwrap();
        let coordinate = locate_track_coordinate_in_range_near(
            track[0],
            &track,
            &cumulative,
            0.0,
            total,
            Some(total),
        )
        .unwrap();

        assert!((coordinate.angle - total).abs() < 1e-12);
    }

    #[test]
    fn plane_is_oriented_monotonic_and_less_than_pi() {
        let track = great_circle();
        let fit = fit_great_circle_plane(&track, &(20..=60).collect::<Vec<_>>()).unwrap();
        assert!(
            fit.horizontal_angles
                .windows(2)
                .all(|pair| pair[1] > pair[0])
        );
        assert!(
            fit.horizontal_angles.last().unwrap() - fit.horizontal_angles.first().unwrap() < PI
        );
        assert!(fit.max_out_of_plane_angle < 1e-12);
    }

    #[test]
    fn plane_rejects_pi_span() {
        let track = (0..=40)
            .map(|index| {
                let angle = index as f64 * PI / 40.0;
                Vec3::new(angle.sin(), 0.0, angle.cos())
            })
            .collect::<Vec<_>>();
        assert!(fit_great_circle_plane(&track, &(0..=40).collect::<Vec<_>>()).is_err());
    }

    #[test]
    fn plane_rejects_a_projected_track_that_folds_back_on_itself() {
        let track = [0.0_f64, 0.2, 0.1, 0.3]
            .into_iter()
            .map(|angle| Vec3::new(angle.sin(), 0.0, angle.cos()))
            .collect::<Vec<_>>();

        assert!(fit_great_circle_plane(&track, &[0, 1, 2, 3]).is_err());
    }

    #[test]
    fn ciddor_matches_standard_dry_air_reference() {
        let result = ciddor_index(288.15, 101_325.0, 0.0, 400.0, 600.0).unwrap();
        assert!((result - 1.000_276_975_860_828_5).abs() < 1e-15);
    }
}
