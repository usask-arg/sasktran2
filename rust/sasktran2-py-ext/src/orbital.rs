use std::collections::{BTreeMap, HashMap, HashSet};

use nalgebra::Vector3;
use ndarray::{Array1, Array2, Array3, s};
use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2, ToPyArray};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

use crate::config::PyConfig;
use crate::geodetic::PyGeodetic;
use crate::geometry::InterpolationMethod as PyInterpolationMethod;
use crate::prelude::IntoPyResult;
use sasktran2_core::orbital as orbital_core;
use sasktran2_core::raytracer::Vec3;
use sasktran2_rs::bindings::atmosphere::Atmosphere;
use sasktran2_rs::bindings::brdf::BrdfKind;
use sasktran2_rs::bindings::config::{Config, ThreadingModel};
use sasktran2_rs::bindings::engine::{
    Engine, LinearizationMode, calculate_prepared_jvps, calculate_prepared_radiances,
    calculate_prepared_vjps, set_2d_refractive_profiles_parallel,
};
use sasktran2_rs::bindings::geodetic::Geodetic;
use sasktran2_rs::bindings::geometry::{Geometry2D, InterpolationMethod};
use sasktran2_rs::bindings::output::{JvpOutput, Output, VjpOutput};
use sasktran2_rs::bindings::prelude::Stokes;
use sasktran2_rs::bindings::viewing_geometry::ViewingGeometry;

type PyTrackProjection<'py> = (
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<bool>>,
);

fn normalize(value: Vector3<f64>, name: &str) -> PyResult<Vector3<f64>> {
    let norm = value.norm();
    if !value.iter().all(|item| item.is_finite()) || !norm.is_finite() || norm == 0.0 {
        return Err(PyValueError::new_err(format!(
            "{name} must contain finite, non-zero vectors"
        )));
    }
    Ok(value / norm)
}

fn row_vector(values: &Array2<f64>, row: usize) -> Vector3<f64> {
    Vector3::new(values[[row, 0]], values[[row, 1]], values[[row, 2]])
}

fn vector_array(value: Vector3<f64>) -> [f64; 3] {
    [value.x, value.y, value.z]
}

fn core_vector(value: Vector3<f64>) -> Vec3 {
    Vec3::new(value.x, value.y, value.z)
}

fn nalgebra_vector(value: Vec3) -> Vector3<f64> {
    Vector3::new(value.x, value.y, value.z)
}

fn interpolation_method(value: &PyInterpolationMethod) -> InterpolationMethod {
    match value {
        PyInterpolationMethod::LinearInterpolation => InterpolationMethod::Linear,
        PyInterpolationMethod::ShellInterpolation => InterpolationMethod::Shell,
        PyInterpolationMethod::LowerInterpolation => InterpolationMethod::Lower,
    }
}

fn time_group_duration_ns(value_s: f64) -> PyResult<i64> {
    let value_ns = value_s * 1e9;
    if !value_s.is_finite()
        || value_s <= 0.0
        || !value_ns.is_finite()
        || value_ns < 0.5
        || value_ns > i64::MAX as f64
    {
        return Err(PyValueError::new_err(
            "time_group_duration_s must be finite, positive, and representable as at least one integer nanosecond",
        ));
    }
    Ok(value_ns.round() as i64)
}

fn geoid_surface_point_from_direction(
    direction: Vector3<f64>,
    equatorial_radius_m: f64,
    flattening_factor: f64,
) -> PyResult<Vector3<f64>> {
    let direction = normalize(direction, "geoid radial direction")?;
    let polar_radius_m = equatorial_radius_m * (1.0 - flattening_factor);
    if !equatorial_radius_m.is_finite()
        || equatorial_radius_m <= 0.0
        || !polar_radius_m.is_finite()
        || polar_radius_m <= 0.0
    {
        return Err(PyValueError::new_err(
            "The geoid equatorial and polar radii must be finite and positive",
        ));
    }
    let inverse_radius_squared = (direction.x * direction.x + direction.y * direction.y)
        / equatorial_radius_m.powi(2)
        + direction.z * direction.z / polar_radius_m.powi(2);
    Ok(direction / inverse_radius_squared.sqrt())
}

#[derive(Clone)]
struct OrbitalGridProvenance {
    vertical_slice_indices: Array1<i64>,
    reference_times_ns: Array1<i64>,
    along_track_angles: Array1<f64>,
    ground_track_ecef_m: Array2<f64>,
    grid_reference_times_ns: Array1<i64>,
}

#[pyclass(unsendable)]
pub struct PyOrbitalPlaneGeometry {
    earth_radius_m: f64,
    surface_radii_m: Array1<f64>,
    altitude_grid_m: Array1<f64>,
    ground_track_ecef_m: Array2<f64>,
    track_directions: Vec<Vector3<f64>>,
    cumulative_angles: Array1<f64>,
    interpolation_method: InterpolationMethod,
    geoid_equatorial_radius_m: Option<f64>,
    geoid_flattening_factor: Option<f64>,
    provenance: Option<OrbitalGridProvenance>,
}

impl PyOrbitalPlaneGeometry {
    fn from_arrays(
        earth_radius_m: f64,
        altitude_grid_m: Array1<f64>,
        mut ground_track_ecef_m: Array2<f64>,
        interpolation_method: InterpolationMethod,
        surface_radii_m: Option<Array1<f64>>,
        geoid_parameters: Option<(f64, f64)>,
        provenance: Option<OrbitalGridProvenance>,
    ) -> PyResult<Self> {
        if !earth_radius_m.is_finite() || earth_radius_m <= 0.0 {
            return Err(PyValueError::new_err(
                "earth_radius_m must be finite and positive",
            ));
        }
        if altitude_grid_m.len() < 2
            || altitude_grid_m.iter().any(|value| !value.is_finite())
            || altitude_grid_m
                .windows(2)
                .into_iter()
                .any(|pair| pair[1] <= pair[0])
        {
            return Err(PyValueError::new_err(
                "altitude_grid_m must contain at least two finite, strictly increasing values",
            ));
        }
        if ground_track_ecef_m.ncols() != 3 || ground_track_ecef_m.nrows() < 2 {
            return Err(PyValueError::new_err(
                "ground_track_ecef_m must have shape (orbital_position, 3) with at least two positions",
            ));
        }
        let mut track_directions = Vec::with_capacity(ground_track_ecef_m.nrows());
        for index in 0..ground_track_ecef_m.nrows() {
            track_directions.push(normalize(
                row_vector(&ground_track_ecef_m, index),
                "ground_track_ecef_m",
            )?);
        }
        let surface_radii_m = if let Some(radii) = surface_radii_m {
            if radii.len() != ground_track_ecef_m.nrows()
                || radii
                    .iter()
                    .any(|radius| !radius.is_finite() || *radius <= 0.0)
            {
                return Err(PyValueError::new_err(
                    "surface_radii_m must contain one finite, positive value per orbital position",
                ));
            }
            radii
        } else {
            Array1::from_elem(ground_track_ecef_m.nrows(), earth_radius_m)
        };
        for (index, direction) in track_directions.iter().enumerate() {
            for component in 0..3 {
                ground_track_ecef_m[[index, component]] =
                    surface_radii_m[index] * direction[component];
            }
        }
        let cumulative_angles = Array1::from_vec(
            orbital_core::cumulative_track_angles(
                &track_directions
                    .iter()
                    .copied()
                    .map(core_vector)
                    .collect::<Vec<_>>(),
            )
            .map_err(PyValueError::new_err)?,
        );

        Ok(Self {
            earth_radius_m: surface_radii_m.sum() / surface_radii_m.len() as f64,
            surface_radii_m,
            altitude_grid_m,
            ground_track_ecef_m,
            track_directions,
            cumulative_angles,
            interpolation_method,
            geoid_equatorial_radius_m: geoid_parameters.map(|value| value.0),
            geoid_flattening_factor: geoid_parameters.map(|value| value.1),
            provenance,
        })
    }
}

#[pymethods]
impl PyOrbitalPlaneGeometry {
    #[new]
    #[pyo3(signature = (
        earth_radius_m,
        altitude_grid_m,
        ground_track_ecef_m,
        interpolation,
        surface_radii_m=None,
    ))]
    fn new(
        earth_radius_m: f64,
        altitude_grid_m: PyReadonlyArray1<f64>,
        ground_track_ecef_m: PyReadonlyArray2<f64>,
        interpolation: PyRef<'_, PyInterpolationMethod>,
        surface_radii_m: Option<PyReadonlyArray1<f64>>,
    ) -> PyResult<Self> {
        let altitude_grid_m = altitude_grid_m.as_array().to_owned();
        let ground_track_ecef_m = ground_track_ecef_m.as_array().to_owned();
        let surface_radii_m = surface_radii_m.map(|radii| radii.as_array().to_owned());
        Self::from_arrays(
            earth_radius_m,
            altitude_grid_m,
            ground_track_ecef_m,
            interpolation_method(&interpolation),
            surface_radii_m,
            None,
            None,
        )
    }

    fn altitudes<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.altitude_grid_m.to_pyarray(py)
    }

    fn ground_track_ecef<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.ground_track_ecef_m.to_pyarray(py)
    }

    fn cumulative_angles<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.cumulative_angles.to_pyarray(py)
    }

    fn surface_radii<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.surface_radii_m.to_pyarray(py)
    }

    fn geodetic_coordinates_degrees<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let Some(equatorial_radius_m) = self.geoid_equatorial_radius_m else {
            return Array2::zeros((0, 2)).to_pyarray(py);
        };
        let flattening_factor = self.geoid_flattening_factor.unwrap();
        let polar_radius_m = equatorial_radius_m * (1.0 - flattening_factor);
        let mut result = Array2::zeros((self.ground_track_ecef_m.nrows(), 2));
        for index in 0..self.ground_track_ecef_m.nrows() {
            let point = row_vector(&self.ground_track_ecef_m, index);
            let horizontal = point.x.hypot(point.y);
            result[[index, 0]] = (point.z * equatorial_radius_m.powi(2))
                .atan2(horizontal * polar_radius_m.powi(2))
                .to_degrees();
            result[[index, 1]] = point.y.atan2(point.x).to_degrees();
        }
        result.to_pyarray(py)
    }

    fn anchor_vertical_slice_indices<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i64>> {
        self.provenance
            .as_ref()
            .map_or_else(
                || Array1::zeros(0),
                |value| value.vertical_slice_indices.clone(),
            )
            .to_pyarray(py)
    }

    fn anchor_reference_times_ns<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i64>> {
        self.provenance
            .as_ref()
            .map_or_else(
                || Array1::zeros(0),
                |value| value.reference_times_ns.clone(),
            )
            .to_pyarray(py)
    }

    fn anchor_along_track_angles<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.provenance
            .as_ref()
            .map_or_else(
                || Array1::zeros(0),
                |value| value.along_track_angles.clone(),
            )
            .to_pyarray(py)
    }

    fn anchor_ground_track_ecef<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.provenance
            .as_ref()
            .map_or_else(
                || Array2::zeros((0, 3)),
                |value| value.ground_track_ecef_m.clone(),
            )
            .to_pyarray(py)
    }

    fn grid_reference_times_ns<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i64>> {
        self.provenance
            .as_ref()
            .map_or_else(
                || Array1::zeros(0),
                |value| value.grid_reference_times_ns.clone(),
            )
            .to_pyarray(py)
    }

    #[pyo3(signature = (ecef_positions_m, order="independent"))]
    fn project_ground_track<'py>(
        &self,
        py: Python<'py>,
        ecef_positions_m: PyReadonlyArray2<f64>,
        order: &str,
    ) -> PyResult<PyTrackProjection<'py>> {
        let positions = ecef_positions_m.as_array().to_owned();
        if positions.ncols() != 3 {
            return Err(PyValueError::new_err(
                "ecef_positions_m must have shape (position, 3)",
            ));
        }
        if !matches!(order, "independent" | "increasing" | "decreasing") {
            return Err(PyValueError::new_err(
                "order must be 'independent', 'increasing', or 'decreasing'",
            ));
        }
        let core_track = self
            .track_directions
            .iter()
            .copied()
            .map(core_vector)
            .collect::<Vec<_>>();
        let cumulative = self.cumulative_angles.as_slice().unwrap();
        let final_angle = *cumulative.last().unwrap();
        let mut minimum_angle = cumulative[0];
        let mut maximum_angle = final_angle;
        let mut segments = Array1::zeros(positions.nrows());
        let mut fractions = Array1::zeros(positions.nrows());
        let mut angles = Array1::zeros(positions.nrows());
        let mut residuals = Array1::zeros(positions.nrows());
        let mut at_track_edge = Array1::from_elem(positions.nrows(), false);
        for index in 0..positions.nrows() {
            let target = normalize(row_vector(&positions, index), "ecef_positions_m")?;
            let coordinate = match order {
                "increasing" => orbital_core::locate_track_coordinate_in_range_near(
                    core_vector(target),
                    &core_track,
                    cumulative,
                    minimum_angle,
                    final_angle,
                    Some(minimum_angle),
                )
                .map_err(PyValueError::new_err)?,
                "decreasing" => orbital_core::locate_track_coordinate_in_range_near(
                    core_vector(target),
                    &core_track,
                    cumulative,
                    cumulative[0],
                    maximum_angle,
                    Some(maximum_angle),
                )
                .map_err(PyValueError::new_err)?,
                _ => orbital_core::locate_track_coordinate(
                    core_vector(target),
                    &core_track,
                    cumulative,
                ),
            };
            if order == "increasing" {
                minimum_angle = coordinate.angle;
            } else if order == "decreasing" {
                maximum_angle = coordinate.angle;
            }
            let projected = interpolate_ground_track(self, coordinate.angle).normalize();
            segments[index] = coordinate.segment as i64;
            fractions[index] = coordinate.fraction;
            angles[index] = coordinate.angle;
            residuals[index] = projected.dot(&target).clamp(-1.0, 1.0).acos();
            at_track_edge[index] = (coordinate.angle - cumulative[0]).abs() <= 1e-12
                || (coordinate.angle - final_angle).abs() <= 1e-12;
        }
        Ok((
            segments.to_pyarray(py),
            fractions.to_pyarray(py),
            angles.to_pyarray(py),
            residuals.to_pyarray(py),
            at_track_edge.to_pyarray(py),
        ))
    }

    fn location_shape(&self) -> (usize, usize) {
        (self.track_directions.len(), self.altitude_grid_m.len())
    }

    #[getter]
    fn earth_radius_m(&self) -> f64 {
        self.earth_radius_m
    }
}

#[pyclass(unsendable)]
pub struct PyOrbitalPlaneViewingGeometry {
    times_ns: Array1<i64>,
    observer_positions_ecef_m: Array2<f64>,
    look_directions_ecef: Array2<f64>,
    vertical_slice_indices: Array1<i64>,
    time_bin_origin_ns: i64,
}

#[derive(Clone)]
struct VerticalSlice {
    index: i64,
    observation_indices: Vec<usize>,
    reference_time_ns: i64,
    first_observation_index: usize,
}

fn ordered_vertical_slices(
    viewing: &PyOrbitalPlaneViewingGeometry,
) -> PyResult<Vec<VerticalSlice>> {
    let mut members = BTreeMap::<i64, Vec<usize>>::new();
    for (observation, &vertical_slice) in viewing.vertical_slice_indices.iter().enumerate() {
        members.entry(vertical_slice).or_default().push(observation);
    }
    let mut slices = members
        .into_iter()
        .map(|(index, observation_indices)| {
            let sum_time = observation_indices
                .iter()
                .map(|&observation| viewing.times_ns[observation] as i128)
                .sum::<i128>();
            VerticalSlice {
                index,
                reference_time_ns: sum_time.div_euclid(observation_indices.len() as i128) as i64,
                first_observation_index: observation_indices[0],
                observation_indices,
            }
        })
        .collect::<Vec<_>>();
    slices.sort_by_key(|slice| (slice.reference_time_ns, slice.first_observation_index));
    if slices.is_empty() {
        return Err(PyValueError::new_err(
            "At least one vertical slice is required",
        ));
    }
    Ok(slices)
}

fn vertical_slice_time_groups(
    viewing: &PyOrbitalPlaneViewingGeometry,
    duration_ns: i64,
) -> PyResult<Vec<Vec<usize>>> {
    if duration_ns <= 0 {
        return Err(PyValueError::new_err(
            "The time-group duration must be positive",
        ));
    }
    let mut grouped = BTreeMap::<i128, Vec<usize>>::new();
    for slice in ordered_vertical_slices(viewing)? {
        let bin = (slice.reference_time_ns as i128 - viewing.time_bin_origin_ns as i128)
            .div_euclid(duration_ns as i128);
        grouped
            .entry(bin)
            .or_default()
            .extend(slice.observation_indices);
    }
    Ok(grouped
        .into_values()
        .map(|mut observations| {
            observations.sort_unstable();
            observations
        })
        .collect())
}

fn interpolate_grid_reference_times(
    grid_angles: &Array1<f64>,
    anchor_angles: &[f64],
    anchor_times_ns: &[i64],
) -> Array1<i64> {
    debug_assert_eq!(anchor_angles.len(), anchor_times_ns.len());
    let mut unique_angles = Vec::new();
    let mut unique_times = Vec::new();
    let mut index = 0;
    while index < anchor_angles.len() {
        let angle = anchor_angles[index];
        let mut sum = anchor_times_ns[index] as i128;
        let mut count = 1_i128;
        index += 1;
        while index < anchor_angles.len() && (anchor_angles[index] - angle).abs() <= 1e-12 {
            sum += anchor_times_ns[index] as i128;
            count += 1;
            index += 1;
        }
        unique_angles.push(angle);
        unique_times.push(sum.div_euclid(count) as i64);
    }
    if unique_angles.len() == 1 {
        return Array1::from_elem(grid_angles.len(), unique_times[0]);
    }
    Array1::from_iter(grid_angles.iter().map(|&angle| {
        if angle <= unique_angles[0] {
            return unique_times[0];
        }
        let final_index = unique_angles.len() - 1;
        if angle >= unique_angles[final_index] {
            return unique_times[final_index];
        }
        let upper = unique_angles.partition_point(|candidate| *candidate < angle);
        let lower = upper - 1;
        let fraction =
            (angle - unique_angles[lower]) / (unique_angles[upper] - unique_angles[lower]);
        let delta = unique_times[upper] as i128 - unique_times[lower] as i128;
        let offset = (delta as f64 * fraction).round() as i128;
        (unique_times[lower] as i128 + offset) as i64
    }))
}

#[pymethods]
impl PyOrbitalPlaneViewingGeometry {
    #[new]
    fn new(
        times_ns: PyReadonlyArray1<i64>,
        observer_positions_ecef_m: PyReadonlyArray2<f64>,
        look_directions_ecef: PyReadonlyArray2<f64>,
        vertical_slice_indices: PyReadonlyArray1<i64>,
        time_bin_origin_ns: i64,
    ) -> PyResult<Self> {
        let times_ns = times_ns.as_array().to_owned();
        let observer_positions_ecef_m = observer_positions_ecef_m.as_array().to_owned();
        let mut look_directions_ecef = look_directions_ecef.as_array().to_owned();
        let vertical_slice_indices = vertical_slice_indices.as_array().to_owned();
        let num_los = times_ns.len();
        if num_los == 0
            || observer_positions_ecef_m.dim() != (num_los, 3)
            || look_directions_ecef.dim() != (num_los, 3)
            || vertical_slice_indices.len() != num_los
        {
            return Err(PyValueError::new_err(
                "times, observer positions, look directions, and vertical slices must have shapes (los,), (los, 3), (los, 3), and (los,)",
            ));
        }
        if times_ns.iter().any(|time| *time == i64::MIN) {
            return Err(PyValueError::new_err("times must not contain NaT"));
        }
        if time_bin_origin_ns == i64::MIN || time_bin_origin_ns > *times_ns.iter().min().unwrap() {
            return Err(PyValueError::new_err(
                "time_bin_origin_ns must be a valid time no later than the earliest observation",
            ));
        }
        if vertical_slice_indices
            .iter()
            .any(|vertical_slice| *vertical_slice < 0)
        {
            return Err(PyValueError::new_err(
                "vertical slice indices must be non-negative",
            ));
        }
        let mut closed = HashSet::new();
        let mut previous = vertical_slice_indices[0];
        for &vertical_slice in vertical_slice_indices.iter().skip(1) {
            if vertical_slice != previous {
                closed.insert(previous);
                if closed.contains(&vertical_slice) {
                    return Err(PyValueError::new_err(
                        "Each vertical slice must occupy one contiguous block of LOS",
                    ));
                }
                previous = vertical_slice;
            }
        }
        for index in 0..num_los {
            normalize(
                row_vector(&observer_positions_ecef_m, index),
                "observer_positions_ecef_m",
            )?;
            let direction = normalize(
                row_vector(&look_directions_ecef, index),
                "look_directions_ecef",
            )?;
            for component in 0..3 {
                look_directions_ecef[[index, component]] = direction[component];
            }
        }
        Ok(Self {
            times_ns,
            observer_positions_ecef_m,
            look_directions_ecef,
            vertical_slice_indices,
            time_bin_origin_ns,
        })
    }

    fn times_ns<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i64>> {
        self.times_ns.to_pyarray(py)
    }

    fn observer_positions_ecef<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.observer_positions_ecef_m.to_pyarray(py)
    }

    fn look_directions_ecef<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.look_directions_ecef.to_pyarray(py)
    }

    fn num_rays(&self) -> usize {
        self.times_ns.len()
    }

    fn construct_atmosphere_geometry(
        &self,
        altitude_grid_m: PyReadonlyArray1<f64>,
        along_track_angle_delta: f64,
        path_padding_angle: f64,
        interpolation: PyRef<'_, PyInterpolationMethod>,
        max_orbital_positions: Option<usize>,
        mut geoid: PyRefMut<'_, PyGeodetic>,
    ) -> PyResult<PyOrbitalPlaneGeometry> {
        let altitude_grid_m = altitude_grid_m.as_array().to_owned();
        if altitude_grid_m.len() < 2
            || altitude_grid_m.iter().any(|value| !value.is_finite())
            || altitude_grid_m
                .windows(2)
                .into_iter()
                .any(|pair| pair[1] <= pair[0])
        {
            return Err(PyValueError::new_err(
                "altitude_grid_m must contain at least two finite, strictly increasing values",
            ));
        }
        if !along_track_angle_delta.is_finite()
            || along_track_angle_delta <= 0.0
            || along_track_angle_delta >= std::f64::consts::PI
        {
            return Err(PyValueError::new_err(
                "along_track_angle_delta must be finite and in (0, pi) radians",
            ));
        }
        if !path_padding_angle.is_finite()
            || !(0.0..std::f64::consts::PI).contains(&path_padding_angle)
        {
            return Err(PyValueError::new_err(
                "path_padding_angle must be finite and in [0, pi) radians",
            ));
        }

        let slices = ordered_vertical_slices(self)?;
        let mut observation_directions = vec![Vector3::zeros(); self.times_ns.len()];
        let mut observation_tangent_altitudes_m = vec![0.0; self.times_ns.len()];
        for (index, observation_direction) in observation_directions.iter_mut().enumerate() {
            let observer = row_vector(&self.observer_positions_ecef_m, index);
            let look = row_vector(&self.look_directions_ecef, index);
            geoid
                .output
                .from_tangent_point(vector_array(observer), vector_array(look))
                .into_pyresult()?;
            observation_tangent_altitudes_m[index] = geoid.output.altitude().into_pyresult()?;
            let latitude = geoid.output.latitude().into_pyresult()?;
            let longitude = geoid.output.longitude().into_pyresult()?;
            geoid
                .output
                .from_lat_lon_alt(latitude, longitude, 0.0)
                .into_pyresult()?;
            let location = geoid.output.location().into_pyresult()?;
            *observation_direction = normalize(
                Vector3::new(location[0], location[1], location[2]),
                "geoid tangent surface location",
            )?;
        }

        let duplicate_tolerance = (along_track_angle_delta * 1e-6).max(1e-10);
        let mut directions: Vec<Vector3<f64>> = Vec::with_capacity(slices.len());
        let mut slice_extents: Vec<f64> = Vec::with_capacity(slices.len());
        let mut slice_tangents: Vec<Vector3<f64>> = Vec::with_capacity(slices.len());
        let mut slice_anchor_directions: Vec<Vector3<f64>> = Vec::with_capacity(slices.len());
        for slice in &slices {
            let direction = normalize(
                slice
                    .observation_indices
                    .iter()
                    .map(|&index| observation_directions[index])
                    .fold(Vector3::zeros(), |sum, value| sum + value),
                "vertical-slice mean tangent location",
            )?;
            slice_anchor_directions.push(direction);
            let extent = slice
                .observation_indices
                .iter()
                .map(|&index| {
                    direction
                        .dot(&observation_directions[index])
                        .clamp(-1.0, 1.0)
                        .acos()
                })
                .fold(0.0_f64, f64::max);
            if extent >= std::f64::consts::FRAC_PI_2 {
                return Err(PyValueError::new_err(format!(
                    "Vertical slice {} spans 90 degrees or more and cannot define one orbital-track anchor",
                    slice.index
                )));
            }
            // A finite-altitude limb image has a systematic horizontal gradient
            // in its tangent locations. Regressing the surface direction against
            // tangent altitude gives that direction without depending on LOS row
            // order. A time gradient handles true vertical scans; only a
            // geometrically degenerate image reaches the later fallbacks.
            let mean_altitude = slice
                .observation_indices
                .iter()
                .map(|&index| observation_tangent_altitudes_m[index])
                .sum::<f64>()
                / slice.observation_indices.len() as f64;
            let altitude_scale = slice
                .observation_indices
                .iter()
                .map(|&index| (observation_tangent_altitudes_m[index] - mean_altitude).powi(2))
                .sum::<f64>()
                .sqrt();
            let mut tangent = if altitude_scale > 0.0 {
                slice
                    .observation_indices
                    .iter()
                    .map(|&index| {
                        (observation_tangent_altitudes_m[index] - mean_altitude)
                            * observation_directions[index]
                    })
                    .fold(Vector3::zeros(), |sum, value| sum + value)
                    / altitude_scale
            } else {
                Vector3::zeros()
            };
            tangent -= tangent.dot(&direction) * direction;
            if !tangent.norm().is_finite() || tangent.norm() < 1e-10 {
                let reference_time = slice.reference_time_ns as i128;
                let time_scale = slice
                    .observation_indices
                    .iter()
                    .map(|&index| (self.times_ns[index] as i128 - reference_time) as f64)
                    .map(|offset| offset * offset)
                    .sum::<f64>()
                    .sqrt();
                tangent = if time_scale > 0.0 {
                    slice
                        .observation_indices
                        .iter()
                        .map(|&index| {
                            ((self.times_ns[index] as i128 - reference_time) as f64)
                                * observation_directions[index]
                        })
                        .fold(Vector3::zeros(), |sum, value| sum + value)
                        / time_scale
                } else {
                    Vector3::zeros()
                };
                tangent -= tangent.dot(&direction) * direction;
            }
            if !tangent.norm().is_finite() || tangent.norm() < 1e-12 {
                tangent = slice
                    .observation_indices
                    .iter()
                    .map(|&index| row_vector(&self.look_directions_ecef, index))
                    .fold(Vector3::zeros(), |sum, value| sum + value);
                tangent -= tangent.dot(&direction) * direction;
            }
            if !tangent.norm().is_finite() || tangent.norm() < 1e-12 {
                let trial = if Vector3::z().dot(&direction).abs() > 0.9 {
                    Vector3::y()
                } else {
                    Vector3::z()
                };
                tangent = trial.cross(&direction);
            }
            tangent = tangent.normalize();
            let is_duplicate = directions.last().is_some_and(|previous: &Vector3<f64>| {
                previous.dot(&direction).clamp(-1.0, 1.0).acos() <= duplicate_tolerance
            });
            if is_duplicate {
                let last = slice_extents.last_mut().unwrap();
                *last = (*last).max(extent);
            } else {
                directions.push(direction);
                slice_extents.push(extent);
                slice_tangents.push(tangent);
            }
        }

        if directions.len() == 1 {
            let radial = directions[0];
            let tangent = slice_tangents[0];
            let padding = (path_padding_angle + slice_extents[0]).max(along_track_angle_delta);
            if padding >= std::f64::consts::PI {
                return Err(PyValueError::new_err(
                    "Vertical-slice extent plus path padding must be less than pi",
                ));
            }
            directions = vec![
                padding.cos() * radial - padding.sin() * tangent,
                radial,
                padding.cos() * radial + padding.sin() * tangent,
            ];
        } else {
            let start_padding = path_padding_angle + slice_extents[0];
            let end_padding = path_padding_angle + slice_extents[slice_extents.len() - 1];
            if start_padding >= std::f64::consts::PI || end_padding >= std::f64::consts::PI {
                return Err(PyValueError::new_err(
                    "Vertical-slice extent plus path padding must be less than pi",
                ));
            }
            let first_length = directions[0].dot(&directions[1]).clamp(-1.0, 1.0).acos();
            let first_tangent =
                (directions[1] - first_length.cos() * directions[0]) / first_length.sin();
            let last_index = directions.len() - 1;
            let last_length = directions[last_index - 1]
                .dot(&directions[last_index])
                .clamp(-1.0, 1.0)
                .acos();
            let last_tangent = (last_length.cos() * directions[last_index]
                - directions[last_index - 1])
                / last_length.sin();
            if start_padding > 0.0 {
                directions.insert(
                    0,
                    start_padding.cos() * directions[0] - start_padding.sin() * first_tangent,
                );
            }
            if end_padding > 0.0 {
                let last = *directions.last().unwrap();
                directions.push(end_padding.cos() * last + end_padding.sin() * last_tangent);
            }
        }

        let cumulative = orbital_core::cumulative_track_angles(
            &directions
                .iter()
                .copied()
                .map(core_vector)
                .collect::<Vec<_>>(),
        )
        .map_err(PyValueError::new_err)?;
        let total_angle = *cumulative.last().unwrap();
        let sampling_delta = along_track_angle_delta;
        let estimated_positions = (total_angle / sampling_delta).ceil() + 1.0;
        if !estimated_positions.is_finite() || estimated_positions > usize::MAX as f64 {
            return Err(PyValueError::new_err(
                "The requested along-track spacing produces an unrepresentable atmosphere grid",
            ));
        }
        if let Some(maximum) = max_orbital_positions
            && estimated_positions > maximum as f64
        {
            return Err(PyValueError::new_err(format!(
                "The requested along-track spacing would generate approximately {} orbital positions, exceeding max_orbital_positions={maximum}",
                estimated_positions as usize
            )));
        }
        let mut sample_angles = Vec::new();
        let mut angle = 0.0;
        while angle < total_angle {
            sample_angles.push(angle);
            angle += sampling_delta;
        }
        if total_angle - sample_angles.last().copied().unwrap_or(0.0) > 1e-12 {
            sample_angles.push(total_angle);
        } else if let Some(last) = sample_angles.last_mut() {
            *last = total_angle;
        }

        let mut sampled_directions = Vec::with_capacity(sample_angles.len());
        let mut segment = 0;
        for angle in sample_angles {
            while segment + 1 < cumulative.len() - 1 && angle > cumulative[segment + 1] {
                segment += 1;
            }
            let segment_length = cumulative[segment + 1] - cumulative[segment];
            let offset = angle - cumulative[segment];
            let tangent = (directions[segment + 1] - segment_length.cos() * directions[segment])
                / segment_length.sin();
            sampled_directions.push(offset.cos() * directions[segment] + offset.sin() * tangent);
        }

        let mut sampled_surface = Array2::zeros((sampled_directions.len(), 3));
        let mut sampled_radii = Array1::zeros(sampled_directions.len());
        for (index, direction) in sampled_directions.into_iter().enumerate() {
            let point = geoid_surface_point_from_direction(
                direction,
                geoid.equatorial_radius_m,
                geoid.flattening_factor,
            )?;
            sampled_radii[index] = point.norm();
            for component in 0..3 {
                sampled_surface[[index, component]] = point[component];
            }
        }
        let mean_radius = sampled_radii.sum() / sampled_radii.len() as f64;
        let mut geometry = PyOrbitalPlaneGeometry::from_arrays(
            mean_radius,
            altitude_grid_m,
            sampled_surface,
            interpolation_method(&interpolation),
            Some(sampled_radii),
            Some((geoid.equatorial_radius_m, geoid.flattening_factor)),
            None,
        )?;

        let core_track = geometry
            .track_directions
            .iter()
            .copied()
            .map(core_vector)
            .collect::<Vec<_>>();
        let cumulative = geometry.cumulative_angles.as_slice().unwrap();
        let mut minimum_angle = cumulative[0];
        let mut anchor_angles = Vec::with_capacity(slices.len());
        let mut anchor_surface = Array2::zeros((slices.len(), 3));
        for (index, direction) in slice_anchor_directions.iter().enumerate() {
            let coordinate = orbital_core::locate_track_coordinate_in_range_near(
                core_vector(*direction),
                &core_track,
                cumulative,
                minimum_angle,
                *cumulative.last().unwrap(),
                Some(minimum_angle),
            )
            .map_err(PyValueError::new_err)?;
            minimum_angle = coordinate.angle;
            anchor_angles.push(coordinate.angle);
            let point = geoid_surface_point_from_direction(
                *direction,
                geoid.equatorial_radius_m,
                geoid.flattening_factor,
            )?;
            for component in 0..3 {
                anchor_surface[[index, component]] = point[component];
            }
        }
        let anchor_times = slices
            .iter()
            .map(|slice| slice.reference_time_ns)
            .collect::<Vec<_>>();
        geometry.provenance = Some(OrbitalGridProvenance {
            vertical_slice_indices: Array1::from_iter(slices.iter().map(|slice| slice.index)),
            reference_times_ns: Array1::from_vec(anchor_times.clone()),
            along_track_angles: Array1::from_vec(anchor_angles.clone()),
            ground_track_ecef_m: anchor_surface,
            grid_reference_times_ns: interpolate_grid_reference_times(
                &geometry.cumulative_angles,
                &anchor_angles,
                &anchor_times,
            ),
        });
        Ok(geometry)
    }
}

#[derive(Clone)]
struct TangentCoordinate {
    segment: usize,
    fraction: f64,
    angle: f64,
}

#[derive(Clone)]
struct OrbitalRayInvariant {
    tangent_coordinate: TangentCoordinate,
    tangent_altitude_m: f64,
    observer_altitude_m: f64,
    surface_radius_m: f64,
    tangent_surface_direction: Vector3<f64>,
    look_direction: Vector3<f64>,
}

#[derive(Clone)]
struct GroupRayPolicy {
    tangent_altitude_m: f64,
    observer_altitude_m: f64,
    horizontal_angle_radians: f64,
    tangent_out_of_plane_angle_radians: f64,
    viewing_azimuth_radians: f64,
    refractive_profile_segment: usize,
    refractive_profile_fraction: f64,
    horizontal_location_residual_radians: f64,
    tangent_direction_error_radians: f64,
}

#[derive(Clone)]
struct GroupLayout {
    observation_indices: Vec<usize>,
    grid_indices: Vec<usize>,
    reference_time_ns: i64,
    reference_position: Vector3<f64>,
    earth_radius_m: f64,
    observation_surface_radii_m: Vec<f64>,
    surface_radius_residuals_m: Vec<f64>,
    maximum_relative_observation_radius_residual: f64,
    maximum_relative_horizontal_scale_residual: f64,
    angular_spacing_relative_residuals: Vec<f64>,
    radius_scale_relative_residuals: Vec<f64>,
    horizontal_distance_relative_residuals: Vec<f64>,
    ray_policies: Vec<GroupRayPolicy>,
    reference_z: Vector3<f64>,
    reference_x: Vector3<f64>,
    horizontal_angles: Vec<f64>,
    max_out_of_plane_angle: f64,
    clipped_start: bool,
    clipped_end: bool,
    padding_angle: f64,
    window_expanded: bool,
}

type PlaneFit = (Vector3<f64>, Vector3<f64>, Vec<f64>, f64);
type GroupReferenceData<'py> = (
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray2<f64>>,
);

#[derive(Clone, PartialEq, Eq)]
struct AtmosphereSignature {
    num_wavel: usize,
    num_location: usize,
    num_legendre_storage: usize,
    volume_mappings: Vec<String>,
    surface_mappings: Vec<String>,
    spatial_surface_rows: usize,
    spatial_surface_derivative: Option<String>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
struct AtmosphereContentKey {
    master_instance_id: u64,
    master_revision: u64,
    surface_generation: u64,
}

#[derive(PartialEq)]
struct LambertianSurfaceState {
    field: Array2<f64>,
    derivative_name: Option<String>,
    spectral_interpolator: Array2<f64>,
    spatial_parameters: bool,
}

const ORBITAL_SURFACE_MAPPING_PREFIX: &str = "__orbital_lambertian_node__";
const ORBITAL_MAPPED_SURFACE_PREFIX: &str = "__orbital_mapped_surface__";

fn local_surface_mapping_name(master_name: &str, local_horizontal: usize) -> String {
    format!("{ORBITAL_SURFACE_MAPPING_PREFIX}{master_name}__{local_horizontal}")
}

fn mapped_surface_output_name(master_name: &str) -> String {
    format!("{ORBITAL_MAPPED_SURFACE_PREFIX}{master_name}")
}

fn observation_ray_invariant(
    observer: Vector3<f64>,
    look: Vector3<f64>,
    geometry: &PyOrbitalPlaneGeometry,
    geoid: Option<&mut Geodetic>,
) -> PyResult<OrbitalRayInvariant> {
    let distance = -observer.dot(&look);
    if !distance.is_finite() || distance <= 0.0 {
        return Err(PyValueError::new_err(
            "OrbitalPlaneViewingGeometry supports limb rays looking toward Earth only",
        ));
    }
    let (tangent_altitude_m, observer_altitude_m, surface_radius_m, surface_direction) =
        if let Some(geoid) = geoid {
            geoid
                .from_tangent_point(vector_array(observer), vector_array(look))
                .into_pyresult()?;
            let tangent_altitude_m = geoid.altitude().into_pyresult()?;
            let latitude = geoid.latitude().into_pyresult()?;
            let longitude = geoid.longitude().into_pyresult()?;
            geoid
                .from_lat_lon_alt(latitude, longitude, 0.0)
                .into_pyresult()?;
            let surface = geoid.location().into_pyresult()?;
            let surface = Vector3::new(surface[0], surface[1], surface[2]);
            let surface_radius_m = surface.norm();
            let surface_direction = normalize(surface, "geoid tangent surface location")?;
            geoid
                .from_xyz(observer.x, observer.y, observer.z)
                .into_pyresult()?;
            let observer_altitude_m = geoid.altitude().into_pyresult()?;
            (
                tangent_altitude_m,
                observer_altitude_m,
                surface_radius_m,
                surface_direction,
            )
        } else {
            let tangent = observer + distance * look;
            let tangent_direction = normalize(tangent, "Cartesian tangent location")?;
            let coordinate = locate_track_coordinate(
                tangent_direction,
                &geometry.track_directions,
                &geometry.cumulative_angles,
            );
            let surface_radius_m = interpolate_surface_radius(geometry, &coordinate);
            (
                tangent.norm() - surface_radius_m,
                observer.norm() - surface_radius_m,
                surface_radius_m,
                tangent_direction,
            )
        };
    if !tangent_altitude_m.is_finite()
        || !observer_altitude_m.is_finite()
        || observer_altitude_m < tangent_altitude_m
    {
        return Err(PyValueError::new_err(
            "The ECEF LOS does not define finite geoid-relative tangent and observer altitudes",
        ));
    }
    let tangent_coordinate = locate_track_coordinate(
        surface_direction,
        &geometry.track_directions,
        &geometry.cumulative_angles,
    );
    Ok(OrbitalRayInvariant {
        tangent_coordinate,
        tangent_altitude_m,
        observer_altitude_m,
        surface_radius_m,
        tangent_surface_direction: surface_direction,
        look_direction: look,
    })
}

fn interpolate_surface_radius(
    geometry: &PyOrbitalPlaneGeometry,
    coordinate: &TangentCoordinate,
) -> f64 {
    let lower = geometry.surface_radii_m[coordinate.segment];
    let upper = geometry.surface_radii_m[coordinate.segment + 1];
    (1.0 - coordinate.fraction) * lower + coordinate.fraction * upper
}

fn interpolate_ground_track(geometry: &PyOrbitalPlaneGeometry, angle: f64) -> Vector3<f64> {
    let upper = geometry
        .cumulative_angles
        .iter()
        .position(|candidate| *candidate >= angle)
        .unwrap_or(geometry.cumulative_angles.len() - 1);
    let segment = upper
        .saturating_sub(1)
        .min(geometry.cumulative_angles.len() - 2);
    let segment_length =
        geometry.cumulative_angles[segment + 1] - geometry.cumulative_angles[segment];
    let offset = (angle - geometry.cumulative_angles[segment]).clamp(0.0, segment_length);
    let fraction = offset / segment_length;
    let tangent = (geometry.track_directions[segment + 1]
        - segment_length.cos() * geometry.track_directions[segment])
        / segment_length.sin();
    let direction = offset.cos() * geometry.track_directions[segment] + offset.sin() * tangent;
    let radius = (1.0 - fraction) * geometry.surface_radii_m[segment]
        + fraction * geometry.surface_radii_m[segment + 1];
    radius * direction
}

fn locate_track_coordinate(
    target: Vector3<f64>,
    track: &[Vector3<f64>],
    cumulative: &Array1<f64>,
) -> TangentCoordinate {
    locate_track_coordinate_near(target, track, cumulative, None)
}

fn locate_track_coordinate_near(
    target: Vector3<f64>,
    track: &[Vector3<f64>],
    cumulative: &Array1<f64>,
    expected_angle: Option<f64>,
) -> TangentCoordinate {
    let coordinate = orbital_core::locate_track_coordinate_near(
        core_vector(target),
        &track.iter().copied().map(core_vector).collect::<Vec<_>>(),
        cumulative.as_slice().unwrap(),
        expected_angle,
    );
    TangentCoordinate {
        segment: coordinate.segment,
        fraction: coordinate.fraction,
        angle: coordinate.angle,
    }
}

fn assign_time_aware_tangent_coordinates(
    geometry: &PyOrbitalPlaneGeometry,
    viewing: &PyOrbitalPlaneViewingGeometry,
    ray_invariants: &mut [OrbitalRayInvariant],
) -> PyResult<()> {
    let slices = ordered_vertical_slices(viewing)?;
    let core_track = geometry
        .track_directions
        .iter()
        .copied()
        .map(core_vector)
        .collect::<Vec<_>>();
    let cumulative = geometry.cumulative_angles.as_slice().unwrap();
    let mut slice_coordinates = HashMap::with_capacity(slices.len());
    let mut minimum_angle = cumulative[0];
    for slice in &slices {
        let anchor = normalize(
            slice
                .observation_indices
                .iter()
                .map(|&index| ray_invariants[index].tangent_surface_direction)
                .fold(Vector3::zeros(), |sum, value| sum + value),
            "vertical-slice mean tangent location",
        )?;
        let coordinate = orbital_core::locate_track_coordinate_in_range(
            core_vector(anchor),
            &core_track,
            cumulative,
            minimum_angle,
            *cumulative.last().unwrap(),
        )
        .map_err(PyValueError::new_err)?;
        minimum_angle = coordinate.angle;
        slice_coordinates.insert(slice.index, coordinate.angle);
    }

    let has_geoid = geometry.geoid_equatorial_radius_m.is_some();
    for slice in slices {
        let expected_angle = slice_coordinates[&slice.index];
        for observation in slice.observation_indices {
            ray_invariants[observation].tangent_coordinate = locate_track_coordinate_near(
                ray_invariants[observation].tangent_surface_direction,
                &geometry.track_directions,
                &geometry.cumulative_angles,
                Some(expected_angle),
            );
            if !has_geoid {
                let coordinate = &ray_invariants[observation].tangent_coordinate;
                let surface_radius_m = interpolate_surface_radius(geometry, coordinate);
                let observer = row_vector(&viewing.observer_positions_ecef_m, observation);
                let look = row_vector(&viewing.look_directions_ecef, observation);
                let distance = -observer.dot(&look);
                let tangent = observer + distance * look;
                ray_invariants[observation].surface_radius_m = surface_radius_m;
                ray_invariants[observation].tangent_altitude_m = tangent.norm() - surface_radius_m;
                ray_invariants[observation].observer_altitude_m =
                    observer.norm() - surface_radius_m;
            }
        }
    }
    Ok(())
}

fn group_ray_policy(
    invariant: &OrbitalRayInvariant,
    grid_indices: &[usize],
    horizontal_angles: &[f64],
    reference_z: Vector3<f64>,
    reference_x: Vector3<f64>,
) -> PyResult<GroupRayPolicy> {
    let reference_y = normalize(
        reference_z.cross(&reference_x),
        "fitted orbital-plane reference-y",
    )?;
    let direction = invariant.tangent_surface_direction;
    let out_of_plane_component = direction.dot(&reference_y).clamp(-1.0, 1.0);
    let projected = direction - out_of_plane_component * reference_y;
    let projected_norm = projected.norm();
    if !projected_norm.is_finite() || projected_norm < 1e-12 {
        return Err(PyValueError::new_err(
            "An orbital tangent location cannot be projected onto its fitted group plane",
        ));
    }
    let projected = projected / projected_norm;
    let mut horizontal_angle_radians = projected
        .dot(&reference_x)
        .atan2(projected.dot(&reference_z));
    let tangent_out_of_plane_angle_radians = (-out_of_plane_component).atan2(projected_norm);

    let local_lower = grid_indices
        .iter()
        .position(|&index| index == invariant.tangent_coordinate.segment)
        .ok_or_else(|| {
            PyRuntimeError::new_err(
                "The selected orbital window does not contain a ray's tangent segment",
            )
        })?;
    if local_lower + 1 >= horizontal_angles.len()
        || grid_indices[local_lower + 1] != invariant.tangent_coordinate.segment + 1
    {
        return Err(PyRuntimeError::new_err(
            "The selected orbital window does not contain both tangent-segment endpoints",
        ));
    }
    let expected_horizontal_angle = (1.0 - invariant.tangent_coordinate.fraction)
        * horizontal_angles[local_lower]
        + invariant.tangent_coordinate.fraction * horizontal_angles[local_lower + 1];
    horizontal_angle_radians +=
        ((expected_horizontal_angle - horizontal_angle_radians) / std::f64::consts::TAU).round()
            * std::f64::consts::TAU;
    if horizontal_angle_radians < horizontal_angles[0] - 1e-12
        || horizontal_angle_radians > horizontal_angles[horizontal_angles.len() - 1] + 1e-12
    {
        return Err(PyRuntimeError::new_err(
            "The transformed tangent location lies outside its selected local 2D window",
        ));
    }
    let local_upper = horizontal_angles
        .partition_point(|angle| *angle <= horizontal_angle_radians)
        .clamp(1, horizontal_angles.len() - 1);
    let local_cell = local_upper - 1;
    if grid_indices[local_cell + 1] != grid_indices[local_cell] + 1 {
        return Err(PyRuntimeError::new_err(
            "The selected orbital window is not contiguous",
        ));
    }
    let refractive_profile_fraction = ((horizontal_angle_radians - horizontal_angles[local_cell])
        / (horizontal_angles[local_cell + 1] - horizontal_angles[local_cell]))
        .clamp(0.0, 1.0);

    // Match Coordinates::unit_vector_from_angles and
    // Coordinates::local_x_y_from_angles. Retaining the out-of-plane angle
    // means the policy reconstructs the geoid surface location rather than
    // flattening the LOS onto the fitted atmospheric plane.
    let in_plane =
        horizontal_angle_radians.cos() * reference_z + horizontal_angle_radians.sin() * reference_x;
    let local_x =
        horizontal_angle_radians.cos() * reference_x - horizontal_angle_radians.sin() * reference_z;
    let reconstructed_direction = tangent_out_of_plane_angle_radians.cos() * in_plane
        - tangent_out_of_plane_angle_radians.sin() * reference_y;
    let local_y = tangent_out_of_plane_angle_radians.cos() * reference_y
        + tangent_out_of_plane_angle_radians.sin() * in_plane;
    let tangent_look = invariant.look_direction
        - invariant.look_direction.dot(&reconstructed_direction) * reconstructed_direction;
    let tangent_look = normalize(
        tangent_look,
        "LOS direction at the transformed tangent point",
    )?;
    let viewing_azimuth_radians = tangent_look.dot(&local_y).atan2(tangent_look.dot(&local_x));
    let tangent_direction_error_radians = reconstructed_direction
        .dot(&direction)
        .clamp(-1.0, 1.0)
        .acos();

    Ok(GroupRayPolicy {
        tangent_altitude_m: invariant.tangent_altitude_m,
        observer_altitude_m: invariant.observer_altitude_m,
        horizontal_angle_radians,
        tangent_out_of_plane_angle_radians,
        viewing_azimuth_radians,
        refractive_profile_segment: grid_indices[local_cell],
        refractive_profile_fraction,
        horizontal_location_residual_radians: horizontal_angle_radians - expected_horizontal_angle,
        tangent_direction_error_radians,
    })
}

fn fit_group_plane(track: &[Vector3<f64>], grid_indices: &[usize]) -> PyResult<PlaneFit> {
    let fit = orbital_core::fit_great_circle_plane(
        &track.iter().copied().map(core_vector).collect::<Vec<_>>(),
        grid_indices,
    )
    .map_err(PyValueError::new_err)?;
    Ok((
        nalgebra_vector(fit.reference_z),
        nalgebra_vector(fit.reference_x),
        fit.horizontal_angles,
        fit.max_out_of_plane_angle,
    ))
}

fn make_layouts(
    geometry: &PyOrbitalPlaneGeometry,
    viewing: &PyOrbitalPlaneViewingGeometry,
    time_group_duration_ns: i64,
    group_padding_angle: f64,
    max_horizontal_scale_residual: Option<f64>,
) -> PyResult<Vec<GroupLayout>> {
    if time_group_duration_ns <= 0 {
        return Err(PyValueError::new_err(
            "time_group_duration_s must be finite and positive",
        ));
    }
    if !group_padding_angle.is_finite()
        || !(0.0..std::f64::consts::PI).contains(&group_padding_angle)
    {
        return Err(PyValueError::new_err(
            "group_padding_angle must be finite and in [0, pi) radians",
        ));
    }
    if max_horizontal_scale_residual.is_some_and(|value| !value.is_finite() || value < 0.0) {
        return Err(PyValueError::new_err(
            "max_horizontal_scale_residual must be finite and non-negative when supplied",
        ));
    }
    let grouped = vertical_slice_time_groups(viewing, time_group_duration_ns)?;

    let mut geoid = match (
        geometry.geoid_equatorial_radius_m,
        geometry.geoid_flattening_factor,
    ) {
        (Some(radius), Some(flattening)) => {
            Some(Geodetic::new(radius, flattening).into_pyresult()?)
        }
        (None, None) => None,
        _ => {
            return Err(PyRuntimeError::new_err(
                "Orbital geometry contains an incomplete geoid definition",
            ));
        }
    };
    let mut ray_invariants = Vec::with_capacity(viewing.times_ns.len());
    for index in 0..viewing.times_ns.len() {
        ray_invariants.push(observation_ray_invariant(
            row_vector(&viewing.observer_positions_ecef_m, index),
            row_vector(&viewing.look_directions_ecef, index),
            geometry,
            geoid.as_mut(),
        )?);
    }
    assign_time_aware_tangent_coordinates(geometry, viewing, &mut ray_invariants)?;
    let tangent_coordinates = ray_invariants
        .iter()
        .map(|ray| ray.tangent_coordinate.clone())
        .collect::<Vec<_>>();

    let mut result = Vec::with_capacity(grouped.len());
    for (group_index, observations) in grouped.into_iter().enumerate() {
        let min_angle = observations
            .iter()
            .map(|&index| tangent_coordinates[index].angle)
            .fold(f64::INFINITY, f64::min)
            - group_padding_angle;
        let max_angle = observations
            .iter()
            .map(|&index| tangent_coordinates[index].angle)
            .fold(f64::NEG_INFINITY, f64::max)
            + group_padding_angle;
        let mut start = geometry
            .cumulative_angles
            .iter()
            .position(|angle| *angle >= min_angle)
            .unwrap_or(geometry.cumulative_angles.len() - 1)
            .saturating_sub(1);
        let mut end = geometry
            .cumulative_angles
            .iter()
            .rposition(|angle| *angle <= max_angle)
            .unwrap_or(0)
            .saturating_add(1)
            .min(geometry.cumulative_angles.len() - 1);
        if end <= start {
            start = start.min(geometry.cumulative_angles.len() - 2);
            end = start + 1;
        }
        let clipped_start = start == 0 && min_angle < geometry.cumulative_angles[0];
        let clipped_end = end == geometry.cumulative_angles.len() - 1
            && max_angle > geometry.cumulative_angles[geometry.cumulative_angles.len() - 1];
        let grid_indices = (start..=end).collect::<Vec<_>>();
        let (reference_z, reference_x, horizontal_angles, max_out_of_plane_angle) =
            fit_group_plane(&geometry.track_directions, &grid_indices)?;
        let sum_time: i128 = observations
            .iter()
            .map(|&index| viewing.times_ns[index] as i128)
            .sum();
        let reference_time_ns = sum_time.div_euclid(observations.len() as i128) as i64;
        let reference_angle = observations
            .iter()
            .map(|&index| tangent_coordinates[index].angle)
            .sum::<f64>()
            / observations.len() as f64;
        let reference_position = interpolate_ground_track(geometry, reference_angle);
        // A local Geometry2D has one spherical radius. Choose the least-squares
        // radius at the actual observation tangent locations, not over the much
        // wider padding window. Tangent altitude is conserved
        // separately by the ray policy; this radius minimizes the remaining
        // angular-to-horizontal-distance scale error over the measured segment.
        let observation_surface_radii_m = observations
            .iter()
            .map(|&index| ray_invariants[index].surface_radius_m)
            .collect::<Vec<_>>();
        let earth_radius_m = observation_surface_radii_m.iter().sum::<f64>()
            / observation_surface_radii_m.len() as f64;
        let surface_radius_residuals_m = observation_surface_radii_m
            .iter()
            .map(|radius| radius - earth_radius_m)
            .collect::<Vec<_>>();
        let maximum_relative_observation_radius_residual = observation_surface_radii_m
            .iter()
            .zip(&surface_radius_residuals_m)
            .map(|(radius, residual)| residual.abs() / radius)
            .fold(0.0_f64, f64::max);
        let mut angular_spacing_relative_residuals =
            Vec::with_capacity(grid_indices.len().saturating_sub(1));
        let mut radius_scale_relative_residuals =
            Vec::with_capacity(grid_indices.len().saturating_sub(1));
        let mut horizontal_distance_relative_residuals =
            Vec::with_capacity(grid_indices.len().saturating_sub(1));
        for local_index in 0..grid_indices.len() - 1 {
            let global_lower = grid_indices[local_index];
            let global_upper = grid_indices[local_index + 1];
            let master_angle =
                geometry.cumulative_angles[global_upper] - geometry.cumulative_angles[global_lower];
            let local_angle = horizontal_angles[local_index + 1] - horizontal_angles[local_index];
            let master_radius = 0.5
                * (geometry.surface_radii_m[global_lower] + geometry.surface_radii_m[global_upper]);
            angular_spacing_relative_residuals.push(local_angle / master_angle - 1.0);
            radius_scale_relative_residuals.push(earth_radius_m / master_radius - 1.0);
            horizontal_distance_relative_residuals
                .push(earth_radius_m * local_angle / (master_radius * master_angle) - 1.0);
        }
        let maximum_relative_horizontal_scale_residual = horizontal_distance_relative_residuals
            .iter()
            .map(|value| value.abs())
            .fold(0.0_f64, f64::max);
        if max_horizontal_scale_residual
            .is_some_and(|limit| maximum_relative_horizontal_scale_residual > limit)
        {
            return Err(PyValueError::new_err(format!(
                "Orbital group {group_index} has maximum horizontal-distance scale residual {maximum_relative_horizontal_scale_residual:.6}, exceeding max_horizontal_scale_residual={:.6}",
                max_horizontal_scale_residual.unwrap()
            )));
        }
        let ray_policies = observations
            .iter()
            .map(|&index| {
                group_ray_policy(
                    &ray_invariants[index],
                    &grid_indices,
                    &horizontal_angles,
                    reference_z,
                    reference_x,
                )
            })
            .collect::<PyResult<Vec<_>>>()?;
        result.push(GroupLayout {
            observation_indices: observations,
            grid_indices,
            reference_time_ns,
            reference_position,
            earth_radius_m,
            observation_surface_radii_m,
            surface_radius_residuals_m,
            maximum_relative_observation_radius_residual,
            maximum_relative_horizontal_scale_residual,
            angular_spacing_relative_residuals,
            radius_scale_relative_residuals,
            horizontal_distance_relative_residuals,
            ray_policies,
            reference_z,
            reference_x,
            horizontal_angles,
            max_out_of_plane_angle,
            clipped_start,
            clipped_end,
            padding_angle: group_padding_angle,
            window_expanded: false,
        });
    }
    Ok(result)
}

/// Return deterministic group reference data before engines are constructed so
/// Python can evaluate its existing solar handler at one time/location per group.
#[pyfunction]
pub fn orbital_group_reference_data<'py>(
    py: Python<'py>,
    geometry: PyRef<'_, PyOrbitalPlaneGeometry>,
    viewing: PyRef<'_, PyOrbitalPlaneViewingGeometry>,
    time_group_duration_s: f64,
    group_padding_angle: f64,
    max_horizontal_scale_residual: Option<f64>,
) -> PyResult<GroupReferenceData<'py>> {
    let duration_ns = time_group_duration_ns(time_group_duration_s)?;
    let groups = make_layouts(
        &geometry,
        &viewing,
        duration_ns,
        group_padding_angle,
        max_horizontal_scale_residual,
    )?;
    let times = Array1::from_iter(groups.iter().map(|group| group.reference_time_ns));
    let mut positions = Array2::zeros((groups.len(), 3));
    let mut geodetic_coordinates = Array2::zeros((groups.len(), 3));
    let mut local_up = Array2::zeros((groups.len(), 3));
    let mut local_north = Array2::zeros((groups.len(), 3));
    let mut local_east = Array2::zeros((groups.len(), 3));
    let mut horizontal_angle_bounds = Array2::zeros((groups.len(), 2));
    let mut reference_geoid = Geodetic::new(
        geometry
            .geoid_equatorial_radius_m
            .unwrap_or(geometry.earth_radius_m),
        geometry.geoid_flattening_factor.unwrap_or(0.0),
    )
    .into_pyresult()?;
    for (index, group) in groups.iter().enumerate() {
        horizontal_angle_bounds[[index, 0]] = group.horizontal_angles[0];
        horizontal_angle_bounds[[index, 1]] =
            group.horizontal_angles[group.horizontal_angles.len() - 1];
        positions.row_mut(index).assign(&Array1::from_vec(
            vector_array(group.reference_position).to_vec(),
        ));
        reference_geoid
            .from_xyz(
                group.reference_position.x,
                group.reference_position.y,
                group.reference_position.z,
            )
            .into_pyresult()?;
        let latitude = reference_geoid.latitude().into_pyresult()?;
        let longitude = reference_geoid.longitude().into_pyresult()?;
        reference_geoid
            .from_lat_lon_alt(latitude, longitude, 0.0)
            .into_pyresult()?;
        geodetic_coordinates[[index, 0]] = latitude;
        geodetic_coordinates[[index, 1]] = longitude;
        geodetic_coordinates[[index, 2]] = 0.0;
        let up = reference_geoid.local_up().into_pyresult()?;
        let south = reference_geoid.local_south().into_pyresult()?;
        let west = reference_geoid.local_west().into_pyresult()?;
        for component in 0..3 {
            local_up[[index, component]] = up[component];
            local_north[[index, component]] = -south[component];
            local_east[[index, component]] = -west[component];
        }
    }
    Ok((
        times.to_pyarray(py),
        positions.to_pyarray(py),
        geodetic_coordinates.to_pyarray(py),
        local_up.to_pyarray(py),
        local_north.to_pyarray(py),
        local_east.to_pyarray(py),
        horizontal_angle_bounds.to_pyarray(py),
    ))
}

struct GroupEngine {
    // Drop the engine before the C++ geometry/viewing objects it references.
    engine: Engine<'static>,
    geometry: Box<Geometry2D>,
    _viewing: Box<ViewingGeometry>,
    layout: GroupLayout,
    local_indices: Vec<usize>,
    atmosphere: Option<Atmosphere>,
    atmosphere_signature: Option<AtmosphereSignature>,
    atmosphere_content_key: Option<AtmosphereContentKey>,
    atmosphere_update_count: usize,
    volume_update_count: usize,
    surface_only_update_count: usize,
    vjp_cotangent: Option<Array3<f64>>,
    last_refractive_profiles: Option<Array2<f64>>,
    geometry_refresh_count: usize,
}

fn build_group_engine(
    config: &'static Config,
    altitude_grid_m: &[f64],
    interpolation_method: InterpolationMethod,
    layout: GroupLayout,
    sun: Vector3<f64>,
) -> PyResult<GroupEngine> {
    let geometry = Box::new(
        Geometry2D::new_ecef(
            layout.earth_radius_m,
            altitude_grid_m.to_vec(),
            layout.horizontal_angles.clone(),
            interpolation_method,
            vector_array(layout.reference_z),
            vector_array(layout.reference_x),
            vector_array(sun),
        )
        .into_pyresult()?,
    );
    let mut viewing = Box::new(ViewingGeometry::new());
    for ray in &layout.ray_policies {
        viewing
            .add_tangent_altitude(
                ray.tangent_altitude_m,
                ray.observer_altitude_m,
                ray.horizontal_angle_radians,
                ray.tangent_out_of_plane_angle_radians,
                ray.viewing_azimuth_radians,
            )
            .into_pyresult()?;
    }
    // The wrappers own stable C++ allocations. The GroupEngine field order
    // ensures Engine is dropped before either referenced wrapper.
    let geometry_ref =
        unsafe { std::mem::transmute::<&Geometry2D, &'static Geometry2D>(geometry.as_ref()) };
    let viewing_ref = unsafe {
        std::mem::transmute::<&ViewingGeometry, &'static ViewingGeometry>(viewing.as_ref())
    };
    let engine = Engine::new_2d(config, geometry_ref, viewing_ref).into_pyresult()?;
    let local_indices = global_location_indices(&layout, altitude_grid_m.len());
    Ok(GroupEngine {
        engine,
        geometry,
        _viewing: viewing,
        layout,
        local_indices,
        atmosphere: None,
        atmosphere_signature: None,
        atmosphere_content_key: None,
        atmosphere_update_count: 0,
        volume_update_count: 0,
        surface_only_update_count: 0,
        vjp_cotangent: None,
        last_refractive_profiles: None,
        geometry_refresh_count: 0,
    })
}

// Every newly built group owns independent geometry, viewing geometry, and
// engine allocations. Moving the group does not move the boxed C++ referents,
// the engine field is dropped before those referents, and no group is accessed
// concurrently during or after this one-time ownership transfer. The C++
// objects have no thread-affine destruction requirement. Mutable atmosphere
// and derivative workspaces are still empty when this wrapper is constructed.
struct SendGroupEngine(GroupEngine);

// SAFETY: the invariants above are established by build_group_engine, and this
// private wrapper is only created for the scoped construction result below.
unsafe impl Send for SendGroupEngine {}

fn build_group_engines_parallel(
    config: &'static Config,
    altitude_grid_m: &[f64],
    interpolation_method: InterpolationMethod,
    layouts: Vec<GroupLayout>,
    sun_vectors: Vec<Vector3<f64>>,
    num_threads: usize,
) -> PyResult<Vec<GroupEngine>> {
    // Construction concurrency is intentionally local. The shared SASKTRAN2
    // Rayon pool must retain the Config thread count used by later radiance and
    // derivative calculations, even when there are fewer groups than threads.
    let thread_pool = sasktran2_rs::util::create_pool(num_threads).into_pyresult()?;

    // Config is immutable for the lifetime of the orbital engine and every
    // child constructor only reads it. Pass its address instead of capturing a
    // non-Send raw-pointer wrapper in the Rayon closure.
    let config_address = config as *const Config as usize;
    let groups = thread_pool.install(|| {
        layouts
            .into_par_iter()
            .zip(sun_vectors.into_par_iter())
            .map(|(layout, sun)| {
                // SAFETY: the constructor's Python Config borrow remains live,
                // is not mutated while the GIL-held constructor is running, and
                // this scoped pool completes before that borrow can end. The
                // completed PyOrbitalPlaneEngine then retains the Config.
                let config = unsafe { &*(config_address as *const Config) };
                build_group_engine(config, altitude_grid_m, interpolation_method, layout, sun)
                    .map(SendGroupEngine)
            })
            .collect::<PyResult<Vec<_>>>()
    })?;
    Ok(groups.into_iter().map(|group| group.0).collect())
}

fn mapping_signature(
    master: &Atmosphere,
    volume_mappings: &[String],
    surface_mappings: &[String],
    lambertian_surface: Option<&LambertianSurfaceState>,
) -> PyResult<AtmosphereSignature> {
    let available_volume = master
        .storage
        .derivative_mapping_names()
        .map_err(PyRuntimeError::new_err)?
        .into_iter()
        .collect::<HashSet<_>>();
    let available_surface = master
        .surface
        .derivative_mapping_names()
        .map_err(PyRuntimeError::new_err)?
        .into_iter()
        .collect::<HashSet<_>>();
    if let Some(name) = volume_mappings
        .iter()
        .find(|name| !available_volume.contains(*name))
    {
        return Err(PyValueError::new_err(format!(
            "Unknown orbital volume derivative mapping {name:?}"
        )));
    }
    if let Some(name) = surface_mappings
        .iter()
        .find(|name| !available_surface.contains(*name))
    {
        return Err(PyValueError::new_err(format!(
            "Unknown orbital surface derivative mapping {name:?}"
        )));
    }
    let mut volume_mappings = volume_mappings.to_vec();
    volume_mappings.sort_unstable();
    volume_mappings.dedup();
    let mut surface_mappings = surface_mappings.to_vec();
    surface_mappings.sort_unstable();
    surface_mappings.dedup();
    let spatial_surface_derivative = lambertian_surface
        .and_then(|surface| surface.derivative_name.as_ref())
        .filter(|name| surface_mappings.contains(name))
        .cloned();
    Ok(AtmosphereSignature {
        num_wavel: master.num_wavel(),
        num_location: master.num_location(),
        num_legendre_storage: master.storage.leg_coeff.dim().0,
        volume_mappings,
        surface_mappings,
        spatial_surface_rows: lambertian_surface.map_or(0, |surface| surface.field.nrows()),
        spatial_surface_derivative,
    })
}

fn all_mapping_names(master: &Atmosphere) -> PyResult<(Vec<String>, Vec<String>)> {
    let volume = master
        .storage
        .derivative_mapping_names()
        .map_err(PyRuntimeError::new_err)?;
    let surface = master
        .surface
        .derivative_mapping_names()
        .map_err(PyRuntimeError::new_err)?;
    Ok((volume, surface))
}

fn global_location_indices(group: &GroupLayout, num_altitudes: usize) -> Vec<usize> {
    group
        .grid_indices
        .iter()
        .flat_map(|horizontal| {
            (0..num_altitudes).map(move |altitude| horizontal * num_altitudes + altitude)
        })
        .collect()
}

fn group_volume_state_matches(group: &GroupEngine, master: &Atmosphere) -> bool {
    let Some(local) = group.atmosphere.as_ref() else {
        return false;
    };
    let Some(signature) = group.atmosphere_signature.as_ref() else {
        return false;
    };
    // Comparing native derivative mappings is substantially more involved
    // than the compact value state. Conservatively rebuild whenever volume
    // mappings are resident; the surface-only retrieval path has none.
    if !signature.volume_mappings.is_empty()
        || !local.storage.solar_irradiance.iter().copied().eq(master
            .storage
            .solar_irradiance
            .iter()
            .copied())
    {
        return false;
    }
    group
        .local_indices
        .iter()
        .enumerate()
        .all(|(local_index, &global_index)| {
            local
                .storage
                .total_extinction
                .row(local_index)
                .iter()
                .copied()
                .eq(master
                    .storage
                    .total_extinction
                    .row(global_index)
                    .iter()
                    .copied())
                && local.storage.ssa.row(local_index).iter().copied().eq(master
                    .storage
                    .ssa
                    .row(global_index)
                    .iter()
                    .copied())
                && local
                    .storage
                    .emission_source
                    .row(local_index)
                    .iter()
                    .copied()
                    .eq(master
                        .storage
                        .emission_source
                        .row(global_index)
                        .iter()
                        .copied())
                && local
                    .storage
                    .leg_coeff
                    .slice(s![.., local_index, ..])
                    .iter()
                    .copied()
                    .eq(master
                        .storage
                        .leg_coeff
                        .slice(s![.., global_index, ..])
                        .iter()
                        .copied())
        })
}

fn update_group_atmosphere(
    group: &mut GroupEngine,
    master: &Atmosphere,
    lambertian_surface: Option<&LambertianSurfaceState>,
    surface_generation: u64,
    volume_mappings: &[String],
    surface_mappings: &[String],
) -> PyResult<()> {
    let signature = mapping_signature(
        master,
        volume_mappings,
        surface_mappings,
        lambertian_surface,
    )?;
    let content_key = AtmosphereContentKey {
        master_instance_id: master.instance_id().into_pyresult()?,
        master_revision: master.revision().into_pyresult()?,
        surface_generation,
    };
    if group.atmosphere_signature.as_ref() == Some(&signature)
        && group.atmosphere_content_key == Some(content_key)
    {
        return Ok(());
    }
    let local_indices = &group.local_indices;
    if group.atmosphere_signature.as_ref() != Some(&signature) {
        let stokes = match master.num_stokes() {
            1 => Stokes::Stokes1,
            3 => Stokes::Stokes3,
            value => {
                return Err(PyValueError::new_err(format!(
                    "Unsupported number of Stokes components: {value}"
                )));
            }
        };
        let num_legendre = signature.num_legendre_storage / stokes.num_legendre();
        group.atmosphere = Some(Atmosphere::new(
            signature.num_wavel,
            local_indices.len(),
            num_legendre,
            !signature.volume_mappings.is_empty() || !signature.surface_mappings.is_empty(),
            true,
            stokes,
            None,
        ));
        group.atmosphere_signature = Some(signature.clone());
        group.atmosphere_content_key = None;
    }
    let volume_changed = !group_volume_state_matches(group, master);
    let local = group.atmosphere.as_mut().unwrap();
    if volume_changed {
        for (local_index, &global_index) in local_indices.iter().enumerate() {
            local
                .storage
                .total_extinction
                .row_mut(local_index)
                .assign(&master.storage.total_extinction.row(global_index));
            local
                .storage
                .ssa
                .row_mut(local_index)
                .assign(&master.storage.ssa.row(global_index));
            local
                .storage
                .emission_source
                .row_mut(local_index)
                .assign(&master.storage.emission_source.row(global_index));
            local
                .storage
                .leg_coeff
                .slice_mut(s![.., local_index, ..])
                .assign(&master.storage.leg_coeff.slice(s![.., global_index, ..]));
        }
        local
            .storage
            .solar_irradiance
            .assign(&master.storage.solar_irradiance);
    }
    if master.surface.brdf_kind != BrdfKind::Lambertian {
        return Err(PyValueError::new_err(
            "OrbitalPlaneEngine v1 supports Lambertian surfaces only",
        ));
    }
    if let Some(surface) = lambertian_surface {
        let mut local_surface =
            Array2::zeros((group.layout.grid_indices.len(), master.num_wavel()));
        local.surface.brdf_args.fill(0.0);
        for (local_horizontal, &global_index) in group.layout.grid_indices.iter().enumerate() {
            local_surface
                .row_mut(local_horizontal)
                .assign(&surface.field.row(global_index));
            local.surface.brdf_args.row_mut(0).scaled_add(
                1.0 / group.layout.grid_indices.len() as f64,
                &surface.field.row(global_index),
            );
        }
        local
            .surface
            .set_spatial_lambertian(local_surface.view())
            .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
    } else {
        local
            .surface
            .clear_spatial_lambertian()
            .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
        local.surface.brdf_args.assign(&master.surface.brdf_args);
    }
    local.surface.emission.assign(&master.surface.emission);

    if volume_changed {
        for name in &signature.volume_mappings {
            let source = master
                .storage
                .get_derivative_mapping(name)
                .map_err(PyRuntimeError::new_err)?;
            let mut destination = local
                .storage
                .get_derivative_mapping(name)
                .map_err(PyRuntimeError::new_err)?;
            for (local_index, &global_index) in local_indices.iter().enumerate() {
                destination
                    .d_extinction()
                    .row_mut(local_index)
                    .assign(&source.d_extinction().row(global_index));
                destination
                    .d_ssa()
                    .row_mut(local_index)
                    .assign(&source.d_ssa().row(global_index));
                destination
                    .d_emission()
                    .row_mut(local_index)
                    .assign(&source.d_emission().row(global_index));
                destination
                    .scat_factor()
                    .row_mut(local_index)
                    .assign(&source.scat_factor().row(global_index));
                destination
                    .d_leg_coeff()
                    .slice_mut(s![.., local_index, ..])
                    .assign(&source.d_leg_coeff().slice(s![.., global_index, ..]));
            }
            destination.set_assign_name(&source.get_assign_name());
            destination.set_interp_dim(&source.get_interp_dim());
            destination.set_log_radiance_space(source.log_radiance_space());
            if let Some(source_interpolator) = source.get_interpolator() {
                let mut gathered =
                    Array2::zeros((local_indices.len(), source_interpolator.ncols()));
                for (local_index, &global_index) in local_indices.iter().enumerate() {
                    gathered
                        .row_mut(local_index)
                        .assign(&source_interpolator.row(global_index));
                }
                destination.set_interpolator(&mut gathered);
            } else {
                // Keep native orbital mappings local. Expanding this identity-like
                // gather to (local locations, all orbital locations) creates a
                // large dense matrix in every resident group engine. The orbital
                // JVP/VJP/Jacobian paths gather or scatter these mappings directly.
                destination.clear_interpolator();
            }
        }
        // Derivative mappings own their phase-function perturbations, while the
        // native source implementations consume a compact scattering-derivative
        // table on AtmosphereGridStorage. The master atmosphere builds that table
        // after registering its constituents, but gathering mappings into this
        // persistent local atmosphere does not do so automatically. Rebuild it
        // after every volume gather so both group membership and copied phase
        // coefficients follow the current atmospheric state.
        local.storage.finalize_scattering_derivatives();
    }
    for name in &signature.surface_mappings {
        if lambertian_surface.and_then(|surface| surface.derivative_name.as_ref()) == Some(name) {
            for local_horizontal in 0..group.layout.grid_indices.len() {
                let local_name = local_surface_mapping_name(name, local_horizontal);
                let mut destination = local
                    .surface
                    .get_derivative_mapping(&local_name)
                    .map_err(PyRuntimeError::new_err)?;
                destination.d_emission().fill(0.0);
                destination.d_brdf().fill(0.0);
                destination.d_brdf().column_mut(local_horizontal).fill(1.0);
                destination.set_interp_dim("wavelength");
                destination.clear_interpolator();
            }
            continue;
        }
        let source = master
            .surface
            .get_derivative_mapping(name)
            .map_err(PyRuntimeError::new_err)?;
        let mut destination = local
            .surface
            .get_derivative_mapping(name)
            .map_err(PyRuntimeError::new_err)?;
        destination.d_emission().assign(&source.d_emission());
        if lambertian_surface.is_some() {
            if source.d_brdf().iter().any(|value| *value != 0.0) {
                return Err(PyValueError::new_err(format!(
                    "Spatial Lambertian surface cannot be combined with BRDF derivative mapping {name:?}"
                )));
            }
            destination.d_brdf().fill(0.0);
        } else {
            destination.d_brdf().assign(&source.d_brdf());
        }
        destination.set_interp_dim(&source.get_interp_dim());
        if let Some(source_interpolator) = source.get_interpolator() {
            let mut interpolator = source_interpolator.to_owned();
            destination.set_interpolator(&mut interpolator);
        }
    }
    if volume_changed {
        local.mark_changed().into_pyresult()?;
        group.volume_update_count += 1;
    } else {
        local.mark_surface_changed().into_pyresult()?;
        group.surface_only_update_count += 1;
    }
    group.atmosphere_content_key = Some(content_key);
    group.atmosphere_update_count += 1;
    Ok(())
}

fn group_vjp_sizes(
    group: &GroupEngine,
    master: &Atmosphere,
    derivative_sizes: &HashMap<String, usize>,
    surface_sizes: &HashMap<String, usize>,
    lambertian_surface: Option<&LambertianSurfaceState>,
    num_orbital_positions: usize,
) -> PyResult<(HashMap<String, usize>, HashMap<String, usize>)> {
    let mut local_derivative_sizes = HashMap::with_capacity(derivative_sizes.len());
    for (name, size) in derivative_sizes {
        let mapping = master
            .storage
            .get_derivative_mapping(name)
            .map_err(PyRuntimeError::new_err)?;
        local_derivative_sizes.insert(
            name.clone(),
            if mapping.get_interpolator().is_some() {
                *size
            } else {
                group.local_indices.len()
            },
        );
    }
    let mut local_surface_sizes = surface_sizes.clone();
    if let Some(surface_state) = lambertian_surface
        && let Some(master_name) = surface_state.derivative_name.as_ref()
        && let Some(size) = surface_sizes.get(master_name)
    {
        let spectral_size = surface_state.spectral_interpolator.ncols();
        let expected = if surface_state.spatial_parameters {
            num_orbital_positions * spectral_size
        } else {
            spectral_size
        };
        if *size != expected {
            return Err(PyValueError::new_err(format!(
                "Orbital surface VJP size for {master_name:?} is {size}; expected {expected}"
            )));
        }
        local_surface_sizes.remove(master_name);
        for local_horizontal in 0..group.layout.grid_indices.len() {
            local_surface_sizes.insert(
                local_surface_mapping_name(master_name, local_horizontal),
                master.num_wavel(),
            );
        }
    }
    Ok((local_derivative_sizes, local_surface_sizes))
}

fn accumulate_group_vjp(
    combined: &mut VjpOutput,
    output: &VjpOutput,
    group: &GroupEngine,
    master: &Atmosphere,
    surface_sizes: &HashMap<String, usize>,
    lambertian_surface: Option<&LambertianSurfaceState>,
) -> PyResult<()> {
    for (local_los, &global_los) in group.layout.observation_indices.iter().enumerate() {
        combined
            .radiance
            .slice_mut(s![.., global_los, ..])
            .assign(&output.radiance.slice(s![.., local_los, ..]));
    }
    for (name, gradient) in &output.derivative_gradients {
        let mapping = master
            .storage
            .get_derivative_mapping(name)
            .map_err(PyRuntimeError::new_err)?;
        let target = combined.derivative_gradients.get_mut(name).unwrap();
        if mapping.get_interpolator().is_some() {
            *target += gradient;
        } else {
            for (local_parameter, &global_parameter) in group.local_indices.iter().enumerate() {
                target[global_parameter] += gradient[local_parameter];
            }
        }
    }
    for (name, gradient) in &output.surface_gradients {
        if name.starts_with(ORBITAL_SURFACE_MAPPING_PREFIX) {
            continue;
        }
        *combined.surface_gradients.get_mut(name).unwrap() += gradient;
    }
    if let Some(surface_state) = lambertian_surface
        && let Some(master_name) = surface_state.derivative_name.as_ref()
        && surface_sizes.contains_key(master_name)
    {
        let target = combined.surface_gradients.get_mut(master_name).unwrap();
        let spectral_size = surface_state.spectral_interpolator.ncols();
        for (local_horizontal, &global_horizontal) in group.layout.grid_indices.iter().enumerate() {
            let local_name = local_surface_mapping_name(master_name, local_horizontal);
            let gradient = output.surface_gradients.get(&local_name).ok_or_else(|| {
                PyRuntimeError::new_err(format!("Missing local spatial surface VJP {local_name:?}"))
            })?;
            for wavelength in 0..master.num_wavel() {
                for spectral_parameter in 0..spectral_size {
                    let parameter = if surface_state.spatial_parameters {
                        global_horizontal * spectral_size + spectral_parameter
                    } else {
                        spectral_parameter
                    };
                    target[parameter] += surface_state.spectral_interpolator
                        [[wavelength, spectral_parameter]]
                        * gradient[wavelength];
                }
            }
        }
    }
    Ok(())
}

fn ciddor_index(
    temperature_k: &Array2<f64>,
    pressure_pa: &Array2<f64>,
    specific_humidity: Option<&Array2<f64>>,
    x_co2: f64,
    wavelength_nm: f64,
) -> PyResult<Array2<f64>> {
    if temperature_k.dim() != pressure_pa.dim() {
        return Err(PyValueError::new_err(
            "pressure_pa and temperature_k must have the same native orbital shape",
        ));
    }
    if specific_humidity.is_some_and(|humidity| humidity.dim() != temperature_k.dim()) {
        return Err(PyValueError::new_err(
            "specific_humidity must have the same native orbital shape as pressure and temperature",
        ));
    }
    let mut result = Array2::zeros(temperature_k.dim());
    for ((index, temperature), pressure) in temperature_k.indexed_iter().zip(pressure_pa.iter()) {
        let humidity = specific_humidity.map_or(0.0, |values| values[index]);
        result[index] =
            orbital_core::ciddor_index(*temperature, *pressure, humidity, x_co2, wavelength_nm)
                .map_err(PyValueError::new_err)?;
    }
    Ok(result)
}

#[pyfunction]
pub fn orbital_ciddor_index<'py>(
    py: Python<'py>,
    temperature_k: PyReadonlyArray2<f64>,
    pressure_pa: PyReadonlyArray2<f64>,
    specific_humidity: Option<PyReadonlyArray2<f64>>,
    x_co2: f64,
    wavelength_nm: f64,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let temperature_k = temperature_k.as_array().to_owned();
    let pressure_pa = pressure_pa.as_array().to_owned();
    let humidity = specific_humidity.map(|value| value.as_array().to_owned());
    Ok(ciddor_index(
        &temperature_k,
        &pressure_pa,
        humidity.as_ref(),
        x_co2,
        wavelength_nm,
    )?
    .to_pyarray(py))
}

#[pyclass(unsendable)]
pub struct PyOrbitalPlaneEngine {
    groups: Vec<GroupEngine>,
    num_observations: usize,
    num_orbital_positions: usize,
    num_altitudes: usize,
    num_stokes: usize,
    refraction_wavelength_nm: f64,
    refraction_co2_ppm: f64,
    los_refraction: bool,
    max_refraction_tangent_altitude_m: f64,
    stream_derivatives: bool,
    lambertian_surface: Option<LambertianSurfaceState>,
    surface_generation: u64,
    state_generation: u64,
    _config: Py<PyConfig>,
    _geometry: Py<PyOrbitalPlaneGeometry>,
    _viewing: Py<PyOrbitalPlaneViewingGeometry>,
}

#[pymethods]
impl PyOrbitalPlaneEngine {
    #[new]
    #[allow(clippy::too_many_arguments)]
    fn new(
        config: PyRef<'_, PyConfig>,
        geometry: PyRef<'_, PyOrbitalPlaneGeometry>,
        viewing: PyRef<'_, PyOrbitalPlaneViewingGeometry>,
        time_group_duration_s: f64,
        group_padding_angle: f64,
        sun_vectors_ecef: PyReadonlyArray2<f64>,
        refraction_wavelength_nm: f64,
        refraction_co2_ppm: f64,
        derivative_execution: &str,
        max_horizontal_scale_residual: Option<f64>,
    ) -> PyResult<Self> {
        let duration_ns = time_group_duration_ns(time_group_duration_s)?;
        if !refraction_wavelength_nm.is_finite() || refraction_wavelength_nm <= 0.0 {
            return Err(PyValueError::new_err(
                "refraction_wavelength_nm must be finite and positive",
            ));
        }
        if !refraction_co2_ppm.is_finite() || refraction_co2_ppm < 0.0 {
            return Err(PyValueError::new_err(
                "refraction_co2_ppm must be finite and non-negative",
            ));
        }
        let stream_derivatives = match derivative_execution {
            "resident" => false,
            "streaming" => true,
            value => {
                return Err(PyValueError::new_err(format!(
                    "derivative_execution must be 'resident' or 'streaming', not {value:?}"
                )));
            }
        };
        let layouts = make_layouts(
            &geometry,
            &viewing,
            duration_ns,
            group_padding_angle,
            max_horizontal_scale_residual,
        )?;
        let sun_vectors = sun_vectors_ecef.as_array().to_owned();
        if sun_vectors.dim() != (layouts.len(), 3) {
            return Err(PyValueError::new_err(format!(
                "sun_vectors_ecef must have shape ({}, 3)",
                layouts.len()
            )));
        }
        if config.config.solar_refraction().into_pyresult()? {
            return Err(PyValueError::new_err(
                "OrbitalPlaneEngine supports line-of-sight refraction only; solar_refraction must be False",
            ));
        }
        let config_ref = unsafe { std::mem::transmute::<&Config, &'static Config>(&config.config) };
        let suns = (0..layouts.len())
            .map(|index| normalize(row_vector(&sun_vectors, index), "sun_vectors_ecef"))
            .collect::<PyResult<Vec<_>>>()?;
        let num_threads = config.config.num_threads().into_pyresult()?;
        let threading_model = config.config.threading_model().into_pyresult()?;
        let groups = if num_threads > 1
            && threading_model == ThreadingModel::Wavelength
            && layouts.len() > 1
        {
            build_group_engines_parallel(
                config_ref,
                geometry.altitude_grid_m.as_slice().ok_or_else(|| {
                    PyRuntimeError::new_err("Orbital altitude grid is not contiguous")
                })?,
                geometry.interpolation_method,
                layouts,
                suns,
                num_threads.min(sun_vectors.nrows()),
            )?
        } else {
            layouts
                .into_iter()
                .zip(suns)
                .map(|(layout, sun)| {
                    build_group_engine(
                        config_ref,
                        geometry.altitude_grid_m.as_slice().ok_or_else(|| {
                            PyRuntimeError::new_err("Orbital altitude grid is not contiguous")
                        })?,
                        geometry.interpolation_method,
                        layout,
                        sun,
                    )
                })
                .collect::<PyResult<Vec<_>>>()?
        };
        let los_refraction = config.config.los_refraction().into_pyresult()?;
        let max_refraction_tangent_altitude_m = config
            .config
            .los_refraction_max_tangent_altitude_m()
            .into_pyresult()?;
        let num_observations = viewing.times_ns.len();
        let num_orbital_positions = geometry.track_directions.len();
        let num_altitudes = geometry.altitude_grid_m.len();
        let num_stokes = config.config.num_stokes().into_pyresult()?;
        Ok(Self {
            groups,
            num_observations,
            num_orbital_positions,
            num_altitudes,
            num_stokes,
            refraction_wavelength_nm,
            refraction_co2_ppm,
            los_refraction,
            max_refraction_tangent_altitude_m,
            stream_derivatives,
            lambertian_surface: None,
            surface_generation: 0,
            state_generation: 0,
            _config: config.into(),
            _geometry: geometry.into(),
            _viewing: viewing.into(),
        })
    }

    fn _prepare_refraction(
        &mut self,
        py: Python<'_>,
        pressure_pa: Option<PyReadonlyArray2<f64>>,
        temperature_k: Option<PyReadonlyArray2<f64>>,
        specific_humidity: Option<PyReadonlyArray2<f64>>,
        refractive_index: Option<PyReadonlyArray2<f64>>,
    ) -> PyResult<()> {
        if !self.los_refraction {
            return Ok(());
        }
        let expected = (self.num_orbital_positions, self.num_altitudes);
        let index = if let Some(index) = refractive_index {
            let index = index.as_array().to_owned();
            if index.dim() != expected
                || index
                    .iter()
                    .any(|value| !value.is_finite() || *value <= 0.0)
            {
                return Err(PyValueError::new_err(format!(
                    "refractive_index must be finite, positive, and have shape {expected:?}"
                )));
            }
            index
        } else {
            let pressure = pressure_pa.ok_or_else(|| {
                PyValueError::new_err(
                    "pressure_pa and temperature_k are required when LOS refraction is enabled without an explicit refractive_index",
                )
            })?;
            let temperature = temperature_k.ok_or_else(|| {
                PyValueError::new_err(
                    "pressure_pa and temperature_k are required when LOS refraction is enabled without an explicit refractive_index",
                )
            })?;
            let pressure = pressure.as_array().to_owned();
            let temperature = temperature.as_array().to_owned();
            if pressure.dim() != expected || temperature.dim() != expected {
                return Err(PyValueError::new_err(format!(
                    "Native orbital pressure and temperature must have shape {expected:?}"
                )));
            }
            let humidity = specific_humidity.map(|value| value.as_array().to_owned());
            if humidity
                .as_ref()
                .is_some_and(|value| value.dim() != expected)
            {
                return Err(PyValueError::new_err(format!(
                    "Native orbital humidity must have shape {expected:?}"
                )));
            }
            ciddor_index(
                &temperature,
                &pressure,
                humidity.as_ref(),
                self.refraction_co2_ppm,
                self.refraction_wavelength_nm,
            )?
        };

        let mut changed_indices = Vec::new();
        let mut changed_profiles = Vec::new();
        for (group_index, group) in self.groups.iter().enumerate() {
            let mut profiles =
                Array2::zeros((group.layout.observation_indices.len(), self.num_altitudes));
            for (local_ray, ray) in group.layout.ray_policies.iter().enumerate() {
                if ray.tangent_altitude_m > self.max_refraction_tangent_altitude_m {
                    profiles.row_mut(local_ray).fill(1.0);
                    continue;
                }
                let left = index.row(ray.refractive_profile_segment);
                let right = index.row(ray.refractive_profile_segment + 1);
                profiles.row_mut(local_ray).assign(
                    &(&left * (1.0 - ray.refractive_profile_fraction)
                        + &right * ray.refractive_profile_fraction),
                );
            }
            if group.last_refractive_profiles.as_ref() != Some(&profiles) {
                changed_indices.push(group_index);
                changed_profiles.push(profiles);
            }
        }

        let composite_threads = self.composite_wavelength_threads(py)?;
        if let Some(num_threads) = composite_threads.filter(|_| changed_indices.len() > 1) {
            let engines = changed_indices
                .iter()
                .map(|&index| &self.groups[index].engine)
                .collect::<Vec<_>>();
            let results =
                set_2d_refractive_profiles_parallel(&engines, &changed_profiles, num_threads)
                    .into_pyresult()?;
            let mut first_error = None;
            for ((group_index, profiles), result) in changed_indices
                .into_iter()
                .zip(changed_profiles)
                .zip(results)
            {
                match result {
                    Ok(()) => {
                        let group = &mut self.groups[group_index];
                        group.last_refractive_profiles = Some(profiles);
                        group.geometry_refresh_count += 1;
                        self.state_generation = self.state_generation.wrapping_add(1);
                    }
                    Err(error) => {
                        if first_error.is_none() {
                            first_error = Some(error);
                        }
                    }
                }
            }
            if let Some(error) = first_error {
                return Err(PyRuntimeError::new_err(error.to_string()));
            }
        } else {
            for (group_index, profiles) in changed_indices.into_iter().zip(changed_profiles) {
                let group = &mut self.groups[group_index];
                group
                    .engine
                    .set_2d_refractive_profiles(profiles.view())
                    .into_pyresult()?;
                group.last_refractive_profiles = Some(profiles);
                group.geometry_refresh_count += 1;
                self.state_generation = self.state_generation.wrapping_add(1);
            }
        }
        Ok(())
    }

    fn _set_lambertian_surface(
        &mut self,
        albedo: Option<PyReadonlyArray2<f64>>,
        derivative_name: Option<String>,
        spectral_interpolator: Option<PyReadonlyArray2<f64>>,
        spatial_parameters: bool,
    ) -> PyResult<()> {
        let Some(albedo) = albedo else {
            if self.lambertian_surface.take().is_some() {
                self.surface_generation = self.surface_generation.wrapping_add(1);
                self.state_generation = self.state_generation.wrapping_add(1);
            }
            return Ok(());
        };
        let albedo = albedo.as_array().to_owned();
        if albedo.nrows() != self.num_orbital_positions
            || albedo.ncols() == 0
            || albedo
                .iter()
                .any(|value| !value.is_finite() || !(0.0..=1.0).contains(value))
        {
            return Err(PyValueError::new_err(format!(
                "Orbital Lambertian albedo must have shape ({}, wavelength), be finite, and lie in [0, 1]",
                self.num_orbital_positions
            )));
        }
        let spectral_interpolator = spectral_interpolator
            .ok_or_else(|| {
                PyValueError::new_err(
                    "Orbital Lambertian albedo requires its spectral interpolator",
                )
            })?
            .as_array()
            .to_owned();
        if spectral_interpolator.nrows() != albedo.ncols()
            || spectral_interpolator.ncols() == 0
            || spectral_interpolator.iter().any(|value| !value.is_finite())
        {
            return Err(PyValueError::new_err(format!(
                "Orbital Lambertian spectral interpolator must have shape ({}, parameter_spectral_point) and be finite",
                albedo.ncols()
            )));
        }
        if !spatial_parameters && spectral_interpolator.ncols() != 1 {
            return Err(PyValueError::new_err(
                "A spatially constant Lambertian parameter must have one spectral parameter",
            ));
        }
        let next_surface = LambertianSurfaceState {
            field: albedo,
            derivative_name,
            spectral_interpolator,
            spatial_parameters,
        };
        if self.lambertian_surface.as_ref() != Some(&next_surface) {
            self.lambertian_surface = Some(next_surface);
            self.surface_generation = self.surface_generation.wrapping_add(1);
            self.state_generation = self.state_generation.wrapping_add(1);
        }
        Ok(())
    }

    fn state_generation(&self) -> u64 {
        self.state_generation
    }

    fn calculate_radiance(
        &mut self,
        py: Python,
        atmosphere: PyRef<'_, crate::atmosphere::PyAtmosphere>,
    ) -> PyResult<Py<crate::output::PyOutput>> {
        self.validate_master(py, &atmosphere.atmosphere)?;
        let master = &atmosphere.atmosphere;
        let (volume_mappings, surface_mappings) = all_mapping_names(master)?;
        let num_stokes = master.num_stokes();
        let mut combined = Output::new(master.num_wavel(), self.num_observations, 0, 0, num_stokes);
        let include_los_optical_depth = self
            ._config
            .borrow(py)
            .config
            .output_los_optical_depth()
            .into_pyresult()?;
        let mut combined_los_optical_depth = include_los_optical_depth
            .then(|| Array2::zeros((master.num_wavel(), self.num_observations)));
        for name in &volume_mappings {
            let size = master
                .storage
                .get_derivative_mapping(name)
                .map_err(PyRuntimeError::new_err)?
                .num_output();
            combined.with_derivative(name, size);
        }
        for name in &surface_mappings {
            if self
                .lambertian_surface
                .as_ref()
                .and_then(|surface| surface.derivative_name.as_ref())
                == Some(name)
            {
                let surface = self.lambertian_surface.as_ref().unwrap();
                let spectral_size = surface.spectral_interpolator.ncols();
                let parameter_size = if surface.spatial_parameters {
                    self.num_orbital_positions * spectral_size
                } else {
                    spectral_size
                };
                combined.with_derivative(&mapped_surface_output_name(name), parameter_size);
            } else {
                combined.with_surface_derivative(name);
            }
        }
        for group in &mut self.groups {
            update_group_atmosphere(
                group,
                master,
                self.lambertian_surface.as_ref(),
                self.surface_generation,
                &volume_mappings,
                &surface_mappings,
            )?;
            let local_indices = &group.local_indices;
            let output = group
                .engine
                .calculate_radiance(group.atmosphere.as_ref().unwrap())
                .into_pyresult()?;
            let local_surface_jacobians = if let Some(surface) = self.lambertian_surface.as_ref() {
                if let Some(master_name) = surface.derivative_name.as_ref() {
                    let num_local_horizontal = group.layout.grid_indices.len();
                    let num_local_los = group.layout.observation_indices.len();
                    let output_size = master.num_wavel() * num_local_los * master.num_stokes();
                    let mut jacobians = Vec::with_capacity(num_local_horizontal);
                    let use_vjp = num_local_horizontal > output_size
                        && group
                            .engine
                            .supports_linearization(LinearizationMode::Vjp)
                            .into_pyresult()?;
                    if !use_vjp {
                        for local_horizontal in 0..num_local_horizontal {
                            let mut tangents = HashMap::with_capacity(1);
                            tangents.insert(
                                local_surface_mapping_name(master_name, local_horizontal),
                                Array1::ones(master.num_wavel()),
                            );
                            jacobians.push(
                                group
                                    .engine
                                    .calculate_jvp(
                                        group.atmosphere.as_ref().unwrap(),
                                        &HashMap::new(),
                                        &tangents,
                                    )
                                    .into_pyresult()?
                                    .jvp
                                    .clone(),
                            );
                        }
                    } else {
                        jacobians.resize_with(num_local_horizontal, || {
                            Array3::zeros((master.num_wavel(), num_local_los, master.num_stokes()))
                        });
                        let mut local_surface_sizes = HashMap::with_capacity(num_local_horizontal);
                        for local_horizontal in 0..num_local_horizontal {
                            local_surface_sizes.insert(
                                local_surface_mapping_name(master_name, local_horizontal),
                                master.num_wavel(),
                            );
                        }
                        let mut cotangent =
                            Array3::zeros((master.num_wavel(), num_local_los, master.num_stokes()));
                        for wavelength in 0..master.num_wavel() {
                            for local_los in 0..num_local_los {
                                for stokes in 0..master.num_stokes() {
                                    cotangent.fill(0.0);
                                    cotangent[[wavelength, local_los, stokes]] = 1.0;
                                    let vjp = group
                                        .engine
                                        .calculate_vjp(
                                            group.atmosphere.as_ref().unwrap(),
                                            &cotangent,
                                            &HashMap::new(),
                                            &local_surface_sizes,
                                        )
                                        .into_pyresult()?;
                                    for (local_horizontal, jacobian) in
                                        jacobians.iter_mut().enumerate()
                                    {
                                        let local_name = local_surface_mapping_name(
                                            master_name,
                                            local_horizontal,
                                        );
                                        let gradient = vjp
                                            .surface_gradients
                                            .get(&local_name)
                                            .ok_or_else(|| {
                                                PyRuntimeError::new_err(format!(
                                                    "Missing local spatial surface VJP {local_name:?}"
                                                ))
                                            })?;
                                        jacobian[[wavelength, local_los, stokes]] =
                                            gradient[wavelength];
                                    }
                                }
                            }
                        }
                    }
                    Some(jacobians)
                } else {
                    None
                }
            } else {
                None
            };
            let local_los_optical_depth =
                include_los_optical_depth.then(|| output.los_optical_depth());
            for (local_los, &global_los) in group.layout.observation_indices.iter().enumerate() {
                combined
                    .radiance
                    .slice_mut(s![.., global_los, ..])
                    .assign(&output.radiance.slice(s![.., local_los, ..]));
                if let (Some(combined_od), Some(local_od)) = (
                    combined_los_optical_depth.as_mut(),
                    local_los_optical_depth.as_ref(),
                ) {
                    combined_od
                        .column_mut(global_los)
                        .assign(&local_od.column(local_los));
                }
                for (name, derivative) in &output.d_radiance {
                    let source = master
                        .storage
                        .get_derivative_mapping(name)
                        .map_err(PyRuntimeError::new_err)?;
                    let target = combined.d_radiance.get_mut(name).unwrap();
                    if source.get_interpolator().is_some() {
                        target
                            .slice_mut(s![.., .., global_los, ..])
                            .assign(&derivative.slice(s![.., .., local_los, ..]));
                    } else {
                        for (local_parameter, &global_parameter) in local_indices.iter().enumerate()
                        {
                            target
                                .slice_mut(s![global_parameter, .., global_los, ..])
                                .scaled_add(
                                    1.0,
                                    &derivative.slice(s![local_parameter, .., local_los, ..]),
                                );
                        }
                    }
                }
                for (name, derivative) in &output.d_radiance_surf {
                    if name.starts_with(ORBITAL_SURFACE_MAPPING_PREFIX) {
                        continue;
                    }
                    combined
                        .d_radiance_surf
                        .get_mut(name)
                        .unwrap()
                        .slice_mut(s![.., global_los, ..])
                        .assign(&derivative.slice(s![.., local_los, ..]));
                }
                if let Some(surface) = self.lambertian_surface.as_ref()
                    && let Some(master_name) = surface.derivative_name.as_ref()
                {
                    let target = combined
                        .d_radiance
                        .get_mut(&mapped_surface_output_name(master_name))
                        .unwrap();
                    let spectral_size = surface.spectral_interpolator.ncols();
                    for local_horizontal in 0..group.layout.grid_indices.len() {
                        let derivative = local_surface_jacobians
                            .as_ref()
                            .and_then(|jacobians| jacobians.get(local_horizontal))
                            .ok_or_else(|| {
                                PyRuntimeError::new_err("Missing local spatial surface JVP block")
                            })?;
                        let global_horizontal = group.layout.grid_indices[local_horizontal];
                        for wavelength in 0..master.num_wavel() {
                            for spectral_parameter in 0..spectral_size {
                                let parameter = if surface.spatial_parameters {
                                    global_horizontal * spectral_size + spectral_parameter
                                } else {
                                    spectral_parameter
                                };
                                let weight =
                                    surface.spectral_interpolator[[wavelength, spectral_parameter]];
                                if weight != 0.0 {
                                    target
                                        .slice_mut(s![parameter, wavelength, global_los, ..])
                                        .scaled_add(
                                            weight,
                                            &derivative.slice(s![wavelength, local_los, ..]),
                                        );
                                }
                            }
                        }
                    }
                }
            }
        }
        if let Some(los_optical_depth) = combined_los_optical_depth {
            combined
                .set_los_optical_depth(los_optical_depth)
                .map_err(PyRuntimeError::new_err)?;
        }
        Py::new(py, crate::output::PyOutput { output: combined })
    }

    fn _calculate_radiance_only(
        &mut self,
        py: Python,
        atmosphere: PyRef<'_, crate::atmosphere::PyAtmosphere>,
    ) -> PyResult<Py<crate::output::PyOutput>> {
        self.validate_master(py, &atmosphere.atmosphere)?;
        let master = &atmosphere.atmosphere;
        let mut combined = Output::new(
            master.num_wavel(),
            self.num_observations,
            0,
            0,
            master.num_stokes(),
        );
        let include_los_optical_depth = self
            ._config
            .borrow(py)
            .config
            .output_los_optical_depth()
            .into_pyresult()?;
        let mut combined_los_optical_depth = include_los_optical_depth
            .then(|| Array2::zeros((master.num_wavel(), self.num_observations)));
        let outputs = self.calculate_group_primals(py, master, &[], &[])?;
        for (group, output) in self.groups.iter().zip(&outputs) {
            let local_los_optical_depth =
                include_los_optical_depth.then(|| output.los_optical_depth());
            for (local_los, &global_los) in group.layout.observation_indices.iter().enumerate() {
                combined
                    .radiance
                    .slice_mut(s![.., global_los, ..])
                    .assign(&output.radiance.slice(s![.., local_los, ..]));
                if let (Some(combined_od), Some(local_od)) = (
                    combined_los_optical_depth.as_mut(),
                    local_los_optical_depth.as_ref(),
                ) {
                    combined_od
                        .column_mut(global_los)
                        .assign(&local_od.column(local_los));
                }
            }
        }
        if let Some(los_optical_depth) = combined_los_optical_depth {
            combined
                .set_los_optical_depth(los_optical_depth)
                .map_err(PyRuntimeError::new_err)?;
        }
        Py::new(py, crate::output::PyOutput { output: combined })
    }

    fn _supports_linearization(&self, mode: u8) -> PyResult<bool> {
        let mode = match mode {
            0 => LinearizationMode::Jacobian,
            1 => LinearizationMode::Jvp,
            2 => LinearizationMode::Vjp,
            _ => {
                return Err(PyValueError::new_err(
                    "linearization mode must be 0, 1, or 2",
                ));
            }
        };
        self.groups
            .iter()
            .map(|group| group.engine.supports_linearization(mode))
            .collect::<Result<Vec<_>, _>>()
            .map(|values| values.into_iter().all(|value| value))
            .into_pyresult()
    }

    fn _linearization_backend(&self, mode: u8) -> PyResult<u8> {
        let mode = match mode {
            0 => LinearizationMode::Jacobian,
            1 => LinearizationMode::Jvp,
            2 => LinearizationMode::Vjp,
            _ => {
                return Err(PyValueError::new_err(
                    "linearization mode must be 0, 1, or 2",
                ));
            }
        };
        let backends = self
            .groups
            .iter()
            .map(|group| group.engine.linearization_backend(mode))
            .collect::<Result<Vec<_>, _>>()
            .into_pyresult()?;
        let backend = backends
            .into_iter()
            .min_by_key(|backend| *backend as i32)
            .ok_or_else(|| PyRuntimeError::new_err("Orbital engine has no groups"))?;
        Ok(backend as u8)
    }

    fn _calculate_jvp(
        &mut self,
        py: Python,
        atmosphere: PyRef<'_, crate::atmosphere::PyAtmosphere>,
        derivative_tangents: &Bound<PyDict>,
        surface_tangents: &Bound<PyDict>,
    ) -> PyResult<Py<crate::output::PyJvpOutput>> {
        self.validate_master(py, &atmosphere.atmosphere)?;
        let volume = array_dict(derivative_tangents)?;
        let surface = array_dict(surface_tangents)?;
        let master = &atmosphere.atmosphere;
        let mut volume_mappings = volume.keys().cloned().collect::<Vec<_>>();
        volume_mappings.sort_unstable();
        let mut surface_mappings = surface.keys().cloned().collect::<Vec<_>>();
        surface_mappings.sort_unstable();
        let mut combined = JvpOutput::new(
            master.num_wavel(),
            self.num_observations,
            master.num_stokes(),
        );
        if volume.is_empty() && surface.is_empty() {
            let outputs = self.calculate_group_primals(py, master, &[], &[])?;
            for (group, output) in self.groups.iter().zip(&outputs) {
                for (local_los, &global_los) in group.layout.observation_indices.iter().enumerate()
                {
                    combined
                        .radiance
                        .slice_mut(s![.., global_los, ..])
                        .assign(&output.radiance.slice(s![.., local_los, ..]));
                }
            }
            return Py::new(py, crate::output::PyJvpOutput { output: combined });
        }
        let composite_threads = self.composite_wavelength_threads(py)?.filter(|_| {
            self.groups.iter().all(|group| {
                matches!(
                    group.engine.linearization_backend(LinearizationMode::Jvp),
                    Ok(sasktran2_rs::bindings::engine::LinearizationBackend::Native)
                )
            })
        });
        let mut outputs = Vec::with_capacity(self.groups.len());
        for group in &mut self.groups {
            update_group_atmosphere(
                group,
                master,
                self.lambertian_surface.as_ref(),
                self.surface_generation,
                &volume_mappings,
                &surface_mappings,
            )?;
            let local_indices = &group.local_indices;
            let mut local_volume = HashMap::with_capacity(volume.len());
            for (name, tangent) in &volume {
                let mapping = master
                    .storage
                    .get_derivative_mapping(name)
                    .map_err(PyRuntimeError::new_err)?;
                if mapping.get_interpolator().is_some() {
                    local_volume.insert(name.clone(), tangent.clone());
                } else {
                    if tangent.len() != master.num_location() {
                        return Err(PyValueError::new_err(format!(
                            "Native orbital tangent {name:?} has length {}; expected {}",
                            tangent.len(),
                            master.num_location()
                        )));
                    }
                    local_volume.insert(
                        name.clone(),
                        Array1::from_iter(local_indices.iter().map(|&index| tangent[index])),
                    );
                }
            }
            let mut local_surface = surface.clone();
            if let Some(surface_state) = self.lambertian_surface.as_ref()
                && let Some(master_name) = surface_state.derivative_name.as_ref()
                && let Some(tangent) = surface.get(master_name)
            {
                local_surface.remove(master_name);
                let spectral_size = surface_state.spectral_interpolator.ncols();
                let expected = if surface_state.spatial_parameters {
                    self.num_orbital_positions * spectral_size
                } else {
                    spectral_size
                };
                if tangent.len() != expected {
                    return Err(PyValueError::new_err(format!(
                        "Orbital surface tangent {master_name:?} has length {}; expected {expected}",
                        tangent.len()
                    )));
                }
                for (local_horizontal, &global_horizontal) in
                    group.layout.grid_indices.iter().enumerate()
                {
                    let mut spectrum = Array1::zeros(master.num_wavel());
                    for wavelength in 0..master.num_wavel() {
                        for spectral_parameter in 0..spectral_size {
                            let parameter = if surface_state.spatial_parameters {
                                global_horizontal * spectral_size + spectral_parameter
                            } else {
                                spectral_parameter
                            };
                            spectrum[wavelength] += surface_state.spectral_interpolator
                                [[wavelength, spectral_parameter]]
                                * tangent[parameter];
                        }
                    }
                    local_surface.insert(
                        local_surface_mapping_name(master_name, local_horizontal),
                        spectrum,
                    );
                }
            }
            outputs.push(if composite_threads.is_some() {
                group
                    .engine
                    .initialize_jvp(
                        group.atmosphere.as_ref().unwrap(),
                        &local_volume,
                        &local_surface,
                    )
                    .into_pyresult()?
            } else {
                group
                    .engine
                    .calculate_jvp(
                        group.atmosphere.as_ref().unwrap(),
                        &local_volume,
                        &local_surface,
                    )
                    .into_pyresult()?
            });
        }
        if let Some(num_threads) = composite_threads {
            let engines = self
                .groups
                .iter()
                .map(|group| &group.engine)
                .collect::<Vec<_>>();
            calculate_prepared_jvps(&engines, &outputs, master.num_wavel(), num_threads)
                .into_pyresult()?;
        }
        for (group, output) in self.groups.iter().zip(&outputs) {
            for (local_los, &global_los) in group.layout.observation_indices.iter().enumerate() {
                combined
                    .radiance
                    .slice_mut(s![.., global_los, ..])
                    .assign(&output.radiance.slice(s![.., local_los, ..]));
                combined
                    .jvp
                    .slice_mut(s![.., global_los, ..])
                    .assign(&output.jvp.slice(s![.., local_los, ..]));
            }
        }
        Py::new(py, crate::output::PyJvpOutput { output: combined })
    }

    fn _calculate_vjp(
        &mut self,
        py: Python,
        atmosphere: PyRef<'_, crate::atmosphere::PyAtmosphere>,
        cotangent: numpy::PyReadonlyArray3<f64>,
        derivative_sizes: &Bound<PyDict>,
        surface_sizes: &Bound<PyDict>,
    ) -> PyResult<Py<crate::output::PyVjpOutput>> {
        self.validate_master(py, &atmosphere.atmosphere)?;
        let cotangent = cotangent.as_array().to_owned();
        let expected = (
            atmosphere.atmosphere.num_wavel(),
            self.num_observations,
            atmosphere.atmosphere.num_stokes(),
        );
        if cotangent.dim() != expected {
            return Err(PyValueError::new_err(format!(
                "Radiance cotangent shape {:?} does not match {expected:?}",
                cotangent.dim()
            )));
        }
        let derivative_sizes = size_dict(derivative_sizes)?;
        let surface_sizes = size_dict(surface_sizes)?;
        let mut volume_mappings = derivative_sizes.keys().cloned().collect::<Vec<_>>();
        volume_mappings.sort_unstable();
        let mut surface_mappings = surface_sizes.keys().cloned().collect::<Vec<_>>();
        surface_mappings.sort_unstable();
        let mut combined = VjpOutput::new(&cotangent);
        for (name, size) in &derivative_sizes {
            combined
                .with_derivative_gradient(name, *size)
                .map_err(PyRuntimeError::new_err)?;
        }
        for (name, size) in &surface_sizes {
            combined
                .with_surface_gradient(name, *size)
                .map_err(PyRuntimeError::new_err)?;
        }
        let master = &atmosphere.atmosphere;
        let composite_threads = if self.stream_derivatives {
            None
        } else {
            self.composite_wavelength_threads(py)?.filter(|_| {
                self.groups.iter().all(|group| {
                    matches!(
                        group.engine.linearization_backend(LinearizationMode::Vjp),
                        Ok(sasktran2_rs::bindings::engine::LinearizationBackend::Native)
                    )
                })
            })
        };
        if let Some(num_threads) = composite_threads {
            let mut outputs = Vec::with_capacity(self.groups.len());
            for group in &mut self.groups {
                update_group_atmosphere(
                    group,
                    master,
                    self.lambertian_surface.as_ref(),
                    self.surface_generation,
                    &volume_mappings,
                    &surface_mappings,
                )?;
                let (local_derivative_sizes, local_surface_sizes) = group_vjp_sizes(
                    group,
                    master,
                    &derivative_sizes,
                    &surface_sizes,
                    self.lambertian_surface.as_ref(),
                    self.num_orbital_positions,
                )?;
                let local_cotangent = group.vjp_cotangent.get_or_insert_with(|| {
                    Array3::zeros((
                        master.num_wavel(),
                        group.layout.observation_indices.len(),
                        master.num_stokes(),
                    ))
                });
                for (local_los, &global_los) in group.layout.observation_indices.iter().enumerate()
                {
                    local_cotangent
                        .slice_mut(s![.., local_los, ..])
                        .assign(&cotangent.slice(s![.., global_los, ..]));
                }
                outputs.push(
                    group
                        .engine
                        .initialize_vjp(
                            group.atmosphere.as_ref().unwrap(),
                            local_cotangent,
                            &local_derivative_sizes,
                            &local_surface_sizes,
                        )
                        .into_pyresult()?,
                );
            }
            let engines = self
                .groups
                .iter()
                .map(|group| &group.engine)
                .collect::<Vec<_>>();
            calculate_prepared_vjps(&engines, &outputs, master.num_wavel(), num_threads)
                .into_pyresult()?;
            for (group, output) in self.groups.iter().zip(&outputs) {
                accumulate_group_vjp(
                    &mut combined,
                    output,
                    group,
                    master,
                    &surface_sizes,
                    self.lambertian_surface.as_ref(),
                )?;
            }
            return Py::new(py, crate::output::PyVjpOutput { output: combined });
        }

        let config = self._config.borrow(py);
        for group in &mut self.groups {
            let streamed_atmosphere = if self.stream_derivatives {
                let resident_atmosphere = group.atmosphere.take();
                let resident_signature = group.atmosphere_signature.take();
                let resident_content_key = group.atmosphere_content_key.take();
                let update_result = update_group_atmosphere(
                    group,
                    master,
                    self.lambertian_surface.as_ref(),
                    self.surface_generation,
                    &volume_mappings,
                    &surface_mappings,
                );
                let derivative_atmosphere = group.atmosphere.take();
                group.atmosphere_signature.take();
                group.atmosphere_content_key.take();
                group.atmosphere = resident_atmosphere;
                group.atmosphere_signature = resident_signature;
                group.atmosphere_content_key = resident_content_key;
                update_result?;
                Some(derivative_atmosphere.ok_or_else(|| {
                    PyRuntimeError::new_err("Failed to build a derivative group atmosphere")
                })?)
            } else {
                update_group_atmosphere(
                    group,
                    master,
                    self.lambertian_surface.as_ref(),
                    self.surface_generation,
                    &volume_mappings,
                    &surface_mappings,
                )?;
                None
            };
            let (local_derivative_sizes, local_surface_sizes) = group_vjp_sizes(
                group,
                master,
                &derivative_sizes,
                &surface_sizes,
                self.lambertian_surface.as_ref(),
                self.num_orbital_positions,
            )?;
            let local_cotangent = group.vjp_cotangent.get_or_insert_with(|| {
                Array3::zeros((
                    master.num_wavel(),
                    group.layout.observation_indices.len(),
                    master.num_stokes(),
                ))
            });
            for (local_los, &global_los) in group.layout.observation_indices.iter().enumerate() {
                local_cotangent
                    .slice_mut(s![.., local_los, ..])
                    .assign(&cotangent.slice(s![.., global_los, ..]));
            }
            let output = if self.stream_derivatives {
                let mut derivative_engine = Engine::new_2d(
                    &config.config,
                    group.geometry.as_ref(),
                    group._viewing.as_ref(),
                )
                .into_pyresult()?;
                if let Some(profiles) = &group.last_refractive_profiles {
                    derivative_engine
                        .set_2d_refractive_profiles(profiles.view())
                        .into_pyresult()?;
                }
                derivative_engine
                    .calculate_vjp(
                        streamed_atmosphere.as_ref().unwrap(),
                        local_cotangent,
                        &local_derivative_sizes,
                        &local_surface_sizes,
                    )
                    .into_pyresult()?
            } else {
                group
                    .engine
                    .calculate_vjp(
                        group.atmosphere.as_ref().unwrap(),
                        local_cotangent,
                        &local_derivative_sizes,
                        &local_surface_sizes,
                    )
                    .into_pyresult()?
            };
            accumulate_group_vjp(
                &mut combined,
                &output,
                group,
                master,
                &surface_sizes,
                self.lambertian_surface.as_ref(),
            )?;
        }
        Py::new(py, crate::output::PyVjpOutput { output: combined })
    }

    fn group_diagnostics(&self, py: Python) -> PyResult<Vec<Py<PyAny>>> {
        let composite_wavelength_threads = self.composite_wavelength_threads(py)?;
        self.groups
            .iter()
            .enumerate()
            .map(|(index, group)| {
                let dict = PyDict::new(py);
                let reference_y = group.layout.reference_z.cross(&group.layout.reference_x);
                let path_edge_usage = group.engine.horizontal_edge_usage().into_pyresult()?;
                let path_edge_clipping_per_los = path_edge_usage
                    .outer_iter()
                    .map(|row| (row[0] != 0, row[1] != 0))
                    .collect::<Vec<_>>();
                let path_edge_clipping = (
                    path_edge_clipping_per_los.iter().any(|value| value.0),
                    path_edge_clipping_per_los.iter().any(|value| value.1),
                );
                dict.set_item("group_index", index)?;
                dict.set_item("observation_indices", &group.layout.observation_indices)?;
                dict.set_item(
                    "vertical_slice_indices",
                    group
                        .layout
                        .observation_indices
                        .iter()
                        .map(|&observation| {
                            self._viewing.borrow(py).vertical_slice_indices[observation]
                        })
                        .collect::<Vec<_>>(),
                )?;
                dict.set_item("grid_indices", &group.layout.grid_indices)?;
                dict.set_item("reference_time_ns", group.layout.reference_time_ns)?;
                dict.set_item(
                    "reference_position_ecef",
                    vector_array(group.layout.reference_position),
                )?;
                dict.set_item("earth_radius_m", group.layout.earth_radius_m)?;
                dict.set_item(
                    "observation_surface_radii_m",
                    &group.layout.observation_surface_radii_m,
                )?;
                dict.set_item(
                    "surface_radius_residuals_m",
                    &group.layout.surface_radius_residuals_m,
                )?;
                dict.set_item(
                    "maximum_absolute_surface_radius_residual_m",
                    group
                        .layout
                        .surface_radius_residuals_m
                        .iter()
                        .map(|value| value.abs())
                        .fold(0.0_f64, f64::max),
                )?;
                dict.set_item(
                    "rms_surface_radius_residual_m",
                    (group
                        .layout
                        .surface_radius_residuals_m
                        .iter()
                        .map(|value| value * value)
                        .sum::<f64>()
                        / group.layout.surface_radius_residuals_m.len() as f64)
                        .sqrt(),
                )?;
                dict.set_item(
                    "maximum_relative_observation_radius_residual",
                    group.layout.maximum_relative_observation_radius_residual,
                )?;
                dict.set_item(
                    "maximum_relative_horizontal_scale_residual",
                    group.layout.maximum_relative_horizontal_scale_residual,
                )?;
                dict.set_item(
                    "angular_spacing_relative_residuals",
                    &group.layout.angular_spacing_relative_residuals,
                )?;
                dict.set_item(
                    "radius_scale_relative_residuals",
                    &group.layout.radius_scale_relative_residuals,
                )?;
                dict.set_item(
                    "horizontal_distance_relative_residuals",
                    &group.layout.horizontal_distance_relative_residuals,
                )?;
                dict.set_item(
                    "rms_relative_horizontal_scale_residual",
                    (group
                        .layout
                        .horizontal_distance_relative_residuals
                        .iter()
                        .map(|value| value * value)
                        .sum::<f64>()
                        / group.layout.horizontal_distance_relative_residuals.len() as f64)
                        .sqrt(),
                )?;
                dict.set_item(
                    "tangent_altitudes_m",
                    group
                        .layout
                        .ray_policies
                        .iter()
                        .map(|ray| ray.tangent_altitude_m)
                        .collect::<Vec<_>>(),
                )?;
                dict.set_item(
                    "observer_altitudes_m",
                    group
                        .layout
                        .ray_policies
                        .iter()
                        .map(|ray| ray.observer_altitude_m)
                        .collect::<Vec<_>>(),
                )?;
                dict.set_item(
                    "ray_horizontal_angles",
                    group
                        .layout
                        .ray_policies
                        .iter()
                        .map(|ray| ray.horizontal_angle_radians)
                        .collect::<Vec<_>>(),
                )?;
                dict.set_item(
                    "ray_out_of_plane_angles",
                    group
                        .layout
                        .ray_policies
                        .iter()
                        .map(|ray| ray.tangent_out_of_plane_angle_radians)
                        .collect::<Vec<_>>(),
                )?;
                dict.set_item(
                    "ray_viewing_azimuths",
                    group
                        .layout
                        .ray_policies
                        .iter()
                        .map(|ray| ray.viewing_azimuth_radians)
                        .collect::<Vec<_>>(),
                )?;
                dict.set_item(
                    "refractive_profile_segments",
                    group
                        .layout
                        .ray_policies
                        .iter()
                        .map(|ray| ray.refractive_profile_segment)
                        .collect::<Vec<_>>(),
                )?;
                dict.set_item(
                    "refractive_profile_fractions",
                    group
                        .layout
                        .ray_policies
                        .iter()
                        .map(|ray| ray.refractive_profile_fraction)
                        .collect::<Vec<_>>(),
                )?;
                dict.set_item(
                    "horizontal_location_residuals_radians",
                    group
                        .layout
                        .ray_policies
                        .iter()
                        .map(|ray| ray.horizontal_location_residual_radians)
                        .collect::<Vec<_>>(),
                )?;
                dict.set_item(
                    "maximum_tangent_direction_error_radians",
                    group
                        .layout
                        .ray_policies
                        .iter()
                        .map(|ray| ray.tangent_direction_error_radians)
                        .fold(0.0_f64, f64::max),
                )?;
                dict.set_item("reference_z", vector_array(group.layout.reference_z))?;
                dict.set_item("reference_x", vector_array(group.layout.reference_x))?;
                dict.set_item("reference_y", vector_array(reference_y))?;
                dict.set_item(
                    "fitted_basis",
                    [
                        vector_array(group.layout.reference_x),
                        vector_array(reference_y),
                        vector_array(group.layout.reference_z),
                    ],
                )?;
                dict.set_item("horizontal_angles", &group.layout.horizontal_angles)?;
                dict.set_item(
                    "horizontal_distances_m",
                    group
                        .layout
                        .horizontal_angles
                        .iter()
                        .map(|angle| angle * group.layout.earth_radius_m)
                        .collect::<Vec<_>>(),
                )?;
                dict.set_item(
                    "max_out_of_plane_angle",
                    group.layout.max_out_of_plane_angle,
                )?;
                dict.set_item(
                    "maximum_out_of_plane_residual",
                    group.layout.max_out_of_plane_angle,
                )?;
                dict.set_item("clipped_start", group.layout.clipped_start)?;
                dict.set_item("clipped_end", group.layout.clipped_end)?;
                dict.set_item("padding_angle", group.layout.padding_angle)?;
                dict.set_item(
                    "edge_clipping",
                    (group.layout.clipped_start, group.layout.clipped_end),
                )?;
                dict.set_item("path_edge_clipping", path_edge_clipping)?;
                dict.set_item("path_edge_clipping_per_los", path_edge_clipping_per_los)?;
                dict.set_item(
                    "uses_extended_horizontal_edge",
                    path_edge_clipping.0 || path_edge_clipping.1,
                )?;
                dict.set_item("window_expanded", group.layout.window_expanded)?;
                dict.set_item("geometry_refresh_count", group.geometry_refresh_count)?;
                dict.set_item(
                    "max_refraction_tangent_altitude_m",
                    self.max_refraction_tangent_altitude_m,
                )?;
                dict.set_item(
                    "refracted_observation_count",
                    if self.los_refraction {
                        group
                            .layout
                            .ray_policies
                            .iter()
                            .filter(|ray| {
                                ray.tangent_altitude_m <= self.max_refraction_tangent_altitude_m
                            })
                            .count()
                    } else {
                        0
                    },
                )?;
                dict.set_item("atmosphere_update_count", group.atmosphere_update_count)?;
                dict.set_item("volume_update_count", group.volume_update_count)?;
                dict.set_item("surface_only_update_count", group.surface_only_update_count)?;
                dict.set_item(
                    "composite_wavelength_scheduler",
                    composite_wavelength_threads.is_some(),
                )?;
                dict.set_item(
                    "composite_wavelength_threads",
                    composite_wavelength_threads.unwrap_or(1),
                )?;
                dict.set_item("engine_identity", group.engine.engine as usize)?;
                dict.set_item(
                    "atmosphere_workspace_identity",
                    group.atmosphere.as_ref().map_or(0, |atmosphere| {
                        atmosphere.storage.total_extinction.as_ptr() as usize
                    }),
                )?;
                dict.set_item(
                    "resident_volume_derivative_mappings",
                    group
                        .atmosphere_signature
                        .as_ref()
                        .map_or_else(Vec::new, |signature| signature.volume_mappings.clone()),
                )?;
                dict.set_item(
                    "resident_surface_derivative_mappings",
                    group
                        .atmosphere_signature
                        .as_ref()
                        .map_or_else(Vec::new, |signature| signature.surface_mappings.clone()),
                )?;
                Ok(dict.into_any().unbind())
            })
            .collect()
    }
}

impl PyOrbitalPlaneEngine {
    fn composite_wavelength_threads(&self, py: Python<'_>) -> PyResult<Option<usize>> {
        let config = self._config.borrow(py);
        let num_threads = config.config.num_threads().into_pyresult()?;
        let threading_model = config.config.threading_model().into_pyresult()?;
        Ok(
            (num_threads > 1 && threading_model == ThreadingModel::Wavelength)
                .then_some(num_threads),
        )
    }

    fn calculate_group_primals(
        &mut self,
        py: Python<'_>,
        master: &Atmosphere,
        volume_mappings: &[String],
        surface_mappings: &[String],
    ) -> PyResult<Vec<Output>> {
        let composite_threads = self.composite_wavelength_threads(py)?;
        for group in &mut self.groups {
            update_group_atmosphere(
                group,
                master,
                self.lambertian_surface.as_ref(),
                self.surface_generation,
                volume_mappings,
                surface_mappings,
            )?;
        }
        if let Some(num_threads) = composite_threads {
            let outputs = self
                .groups
                .iter()
                .map(|group| {
                    group
                        .engine
                        .initialize_radiance_only(group.atmosphere.as_ref().unwrap())
                        .into_pyresult()
                })
                .collect::<PyResult<Vec<_>>>()?;
            let engines = self
                .groups
                .iter()
                .map(|group| &group.engine)
                .collect::<Vec<_>>();
            calculate_prepared_radiances(&engines, &outputs, master.num_wavel(), num_threads)
                .into_pyresult()?;
            Ok(outputs)
        } else {
            self.groups
                .iter()
                .map(|group| {
                    group
                        .engine
                        .calculate_radiance_with_mappings(
                            group.atmosphere.as_ref().unwrap(),
                            &[],
                            &[],
                        )
                        .into_pyresult()
                })
                .collect()
        }
    }

    fn validate_master(&self, py: Python<'_>, master: &Atmosphere) -> PyResult<()> {
        let current_num_stokes = self
            ._config
            .borrow(py)
            .config
            .num_stokes()
            .into_pyresult()?;
        if current_num_stokes != self.num_stokes {
            return Err(PyValueError::new_err(
                "The Config num_stokes changed after OrbitalPlaneEngine construction; construct a new engine",
            ));
        }
        if master.num_stokes() != self.num_stokes {
            return Err(PyValueError::new_err(format!(
                "Atmosphere has {} Stokes components; OrbitalPlaneEngine requires {}",
                master.num_stokes(),
                self.num_stokes
            )));
        }
        if master.num_location() != self.num_orbital_positions * self.num_altitudes {
            return Err(PyValueError::new_err(format!(
                "Atmosphere has {} locations; OrbitalPlaneGeometry requires {}",
                master.num_location(),
                self.num_orbital_positions * self.num_altitudes
            )));
        }
        if self.los_refraction
            && self
                .groups
                .iter()
                .any(|group| group.last_refractive_profiles.is_none())
        {
            return Err(PyValueError::new_err(
                "Prepare orbital refractive profiles before calculating radiance",
            ));
        }
        Ok(())
    }
}

fn array_dict(dict: &Bound<PyDict>) -> PyResult<HashMap<String, Array1<f64>>> {
    let mut result = HashMap::new();
    for (name, value) in dict.iter() {
        let name = name.extract::<String>()?;
        let value = value.extract::<PyReadonlyArray1<f64>>()?;
        result.insert(name, value.as_array().to_owned());
    }
    Ok(result)
}

fn size_dict(dict: &Bound<PyDict>) -> PyResult<HashMap<String, usize>> {
    let mut result = HashMap::new();
    for (name, value) in dict.iter() {
        result.insert(name.extract::<String>()?, value.extract::<usize>()?);
    }
    Ok(result)
}
