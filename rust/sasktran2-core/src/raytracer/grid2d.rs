//! Structured spherical grid varying in altitude and one horizontal angle.

use super::grid::{
    BoundaryTag, CellId, GridError, InterpolationMethod, InterpolationStencil, MediumGrid,
};
use super::primitive::{Intersection, Primitive, PrimitiveId};
use super::vec3::Vec3;

/// Orthonormal reference basis used to measure the horizontal angle.
///
/// Angles rotate from `reference_z` toward `reference_x`; `reference_y` is the
/// common axis contained by every angular boundary half-plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AngularBasis {
    reference_x: Vec3,
    reference_y: Vec3,
    reference_z: Vec3,
}

impl AngularBasis {
    pub const CANONICAL: Self = Self {
        reference_x: Vec3::UNIT_X,
        reference_y: Vec3::UNIT_Y,
        reference_z: Vec3::UNIT_Z,
    };

    pub fn new(reference_z: Vec3, reference_x: Vec3) -> Result<Self, GridError> {
        const BASIS_EPSILON: f64 = 1e-12;

        if !reference_z.is_finite()
            || !reference_x.is_finite()
            || reference_z.norm() <= BASIS_EPSILON
            || reference_x.norm() <= BASIS_EPSILON
        {
            return Err(GridError::InvalidAngularBasis);
        }

        let reference_z = reference_z.normalized();
        if !reference_z.is_finite() || reference_z.norm() <= BASIS_EPSILON {
            return Err(GridError::InvalidAngularBasis);
        }
        let reference_x = reference_x - reference_z * reference_x.dot(reference_z);
        if !reference_x.is_finite() || reference_x.norm() <= BASIS_EPSILON {
            return Err(GridError::InvalidAngularBasis);
        }
        let reference_x = reference_x.normalized();
        let reference_y = reference_z.cross(reference_x).normalized();
        if !reference_x.is_finite()
            || !reference_y.is_finite()
            || reference_y.norm() <= BASIS_EPSILON
        {
            return Err(GridError::InvalidAngularBasis);
        }

        Ok(Self {
            reference_x,
            reference_y,
            reference_z,
        })
    }

    #[inline(always)]
    pub fn reference_x(self) -> Vec3 {
        self.reference_x
    }

    #[inline(always)]
    pub fn reference_y(self) -> Vec3 {
        self.reference_y
    }

    #[inline(always)]
    pub fn reference_z(self) -> Vec3 {
        self.reference_z
    }

    #[inline(always)]
    pub fn radial_direction(self, angle: f64) -> Vec3 {
        angle.sin() * self.reference_x + angle.cos() * self.reference_z
    }

    #[inline(always)]
    pub fn angular_normal(self, angle: f64) -> Vec3 {
        angle.cos() * self.reference_x - angle.sin() * self.reference_z
    }

    #[inline(always)]
    fn raw_angle(self, point: Vec3) -> f64 {
        point
            .dot(self.reference_x)
            .atan2(point.dot(self.reference_z))
    }
}

/// Finite, non-periodic spherical grid in altitude and one horizontal angle.
///
/// The angular span is restricted to less than pi radians so the domain is a
/// single unambiguous wedge. Flattened location storage is altitude-fastest:
/// `location_index = horizontal_index * num_altitudes + altitude_index`.
/// The horizontal coordinate is singular on the `reference_y` axis; a field
/// used on that axis must therefore be independent of horizontal angle.
#[derive(Debug, Clone, PartialEq)]
pub struct StructuredGrid2D {
    earth_radius: f64,
    altitudes: Vec<f64>,
    horizontal_angles: Vec<f64>,
    altitude_interpolation: InterpolationMethod,
    basis: AngularBasis,
}

impl StructuredGrid2D {
    pub fn new(
        earth_radius: f64,
        altitudes: Vec<f64>,
        horizontal_angles: Vec<f64>,
        altitude_interpolation: InterpolationMethod,
    ) -> Result<Self, GridError> {
        Self::new_with_basis(
            earth_radius,
            altitudes,
            horizontal_angles,
            altitude_interpolation,
            AngularBasis::CANONICAL,
        )
    }

    pub fn new_with_basis(
        earth_radius: f64,
        altitudes: Vec<f64>,
        horizontal_angles: Vec<f64>,
        altitude_interpolation: InterpolationMethod,
        basis: AngularBasis,
    ) -> Result<Self, GridError> {
        if !earth_radius.is_finite() || earth_radius < 0.0 {
            return Err(GridError::InvalidEarthRadius);
        }
        validate_altitudes(&altitudes)?;
        if let Some((index, _)) = altitudes.iter().enumerate().find(|(_, altitude)| {
            let radius = earth_radius + **altitude;
            !radius.is_finite() || radius <= 0.0
        }) {
            return Err(GridError::InvalidRadius { index });
        }
        validate_horizontal_angles(&horizontal_angles)?;

        Ok(Self {
            earth_radius,
            altitudes,
            horizontal_angles,
            altitude_interpolation,
            basis,
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
    pub fn horizontal_angles(&self) -> &[f64] {
        &self.horizontal_angles
    }

    #[inline(always)]
    pub fn altitude_interpolation(&self) -> InterpolationMethod {
        self.altitude_interpolation
    }

    #[inline(always)]
    pub fn basis(&self) -> AngularBasis {
        self.basis
    }

    #[inline(always)]
    pub fn num_locations(&self) -> usize {
        self.altitudes.len() * self.horizontal_angles.len()
    }

    #[inline(always)]
    pub fn location_shape(&self) -> (usize, usize) {
        (self.horizontal_angles.len(), self.altitudes.len())
    }

    #[inline(always)]
    pub fn num_cells(&self) -> usize {
        (self.horizontal_angles.len() - 1) * (self.altitudes.len() - 1)
    }

    #[inline(always)]
    pub fn location_index(&self, altitude_index: usize, horizontal_index: usize) -> usize {
        horizontal_index * self.altitudes.len() + altitude_index
    }

    #[inline(always)]
    pub fn altitude_at(&self, point: Vec3) -> f64 {
        point.norm() - self.earth_radius
    }

    pub fn horizontal_angle_at(&self, point: Vec3) -> f64 {
        let raw = self.basis.raw_angle(point);
        let center = (self.horizontal_angles[0]
            + self.horizontal_angles[self.horizontal_angles.len() - 1])
            / 2.0;
        unwrap_angle_near(raw, center)
    }

    pub fn locate_indices(&self, point: Vec3, epsilon: f64) -> Option<(usize, usize)> {
        let altitude_index = locate_axis_cell(&self.altitudes, self.altitude_at(point), epsilon)?;
        let horizontal_index = locate_axis_cell(
            &self.horizontal_angles,
            self.horizontal_angle_at(point),
            epsilon,
        )?;
        Some((altitude_index, horizontal_index))
    }

    pub fn interpolation_weights_at(&self, point: Vec3) -> InterpolationStencil {
        const LOCATION_EPSILON: f64 = 1e-8;

        let altitude = self.altitude_at(point);
        let horizontal_angle = self.horizontal_angle_at(point);
        if locate_axis_cell(&self.altitudes, altitude, LOCATION_EPSILON).is_none()
            || locate_axis_cell(&self.horizontal_angles, horizontal_angle, LOCATION_EPSILON)
                .is_none()
        {
            return InterpolationStencil::new();
        }

        let altitude_weights = axis_interpolation_weights(
            &self.altitudes,
            altitude,
            self.altitude_interpolation,
            LOCATION_EPSILON,
        );
        let horizontal_weights = axis_interpolation_weights(
            &self.horizontal_angles,
            horizontal_angle,
            InterpolationMethod::Linear,
            LOCATION_EPSILON,
        );

        let mut result = InterpolationStencil::new();
        for horizontal in horizontal_weights.iter() {
            for altitude in altitude_weights.iter() {
                result.push(
                    self.location_index(altitude.index, horizontal.index),
                    altitude.weight * horizontal.weight,
                );
            }
        }
        result
    }

    pub fn primitives(&self) -> Vec<Primitive> {
        let mut primitives =
            Vec::with_capacity(self.altitudes.len() + self.horizontal_angles.len());
        for (index, altitude) in self.altitudes.iter().enumerate() {
            let boundary = if index == 0 {
                BoundaryTag::Surface {
                    altitude_index: index,
                }
            } else if index == self.altitudes.len() - 1 {
                BoundaryTag::TopOfAtmosphere {
                    altitude_index: index,
                }
            } else {
                BoundaryTag::Altitude { index }
            };
            primitives.push(Primitive::sphere(
                PrimitiveId(index),
                Vec3::ZERO,
                self.earth_radius + altitude,
                boundary,
            ));
        }
        for (index, angle) in self.horizontal_angles.iter().enumerate() {
            primitives.push(Primitive::angular_plane(
                PrimitiveId(self.altitudes.len() + index),
                self.basis.angular_normal(*angle),
                self.basis.radial_direction(*angle),
                BoundaryTag::Horizontal { index },
            ));
        }
        primitives
    }
}

impl MediumGrid for StructuredGrid2D {
    fn locate_cell(&self, point: Vec3) -> Option<CellId> {
        self.locate_indices(point, 1e-8)
            .map(|(altitude_index, horizontal_index)| CellId::Structured2D {
                altitude_index,
                horizontal_index,
            })
    }

    fn interpolation_weights(&self, point: Vec3) -> InterpolationStencil {
        self.interpolation_weights_at(point)
    }

    fn classify_intersection(&self, intersection: &Intersection) -> BoundaryTag {
        intersection.boundary
    }
}

fn validate_altitudes(altitudes: &[f64]) -> Result<(), GridError> {
    if altitudes.len() < 2 {
        return Err(GridError::TooFewAltitudes);
    }
    for (index, altitude) in altitudes.iter().enumerate() {
        if !altitude.is_finite() {
            return Err(GridError::NonFiniteAltitude { index });
        }
    }
    for index in 1..altitudes.len() {
        if altitudes[index] <= altitudes[index - 1] {
            return Err(GridError::NonIncreasingAltitude { index });
        }
    }
    Ok(())
}

fn validate_horizontal_angles(horizontal_angles: &[f64]) -> Result<(), GridError> {
    if horizontal_angles.len() < 2 {
        return Err(GridError::TooFewHorizontalAngles);
    }
    for (index, angle) in horizontal_angles.iter().enumerate() {
        if !angle.is_finite() {
            return Err(GridError::NonFiniteHorizontalAngle { index });
        }
    }
    for index in 1..horizontal_angles.len() {
        if horizontal_angles[index] <= horizontal_angles[index - 1] {
            return Err(GridError::NonIncreasingHorizontalAngle { index });
        }
    }
    if horizontal_angles[horizontal_angles.len() - 1] - horizontal_angles[0] >= std::f64::consts::PI
    {
        return Err(GridError::HorizontalSpanTooWide);
    }
    Ok(())
}

fn locate_axis_cell(values: &[f64], value: f64, epsilon: f64) -> Option<usize> {
    if value < values[0] - epsilon || value > values[values.len() - 1] + epsilon {
        return None;
    }
    if let Some(index) = values
        .iter()
        .position(|boundary| (value - boundary).abs() <= epsilon)
    {
        return Some(index.min(values.len() - 2));
    }
    let upper = values.partition_point(|boundary| *boundary <= value);
    if upper == 0 || upper >= values.len() {
        None
    } else {
        Some(upper - 1)
    }
}

fn axis_interpolation_weights(
    values: &[f64],
    value: f64,
    interpolation: InterpolationMethod,
    epsilon: f64,
) -> InterpolationStencil {
    let mut stencil = InterpolationStencil::new();
    if let Some(index) = values
        .iter()
        .position(|boundary| (value - boundary).abs() <= epsilon)
    {
        stencil.push(index, 1.0);
        return stencil;
    }
    if value <= values[0] {
        stencil.push(0, 1.0);
        return stencil;
    }
    let last = values.len() - 1;
    if value >= values[last] {
        stencil.push(last, 1.0);
        return stencil;
    }

    let upper = values.partition_point(|boundary| *boundary <= value);
    let lower = upper - 1;
    match interpolation {
        InterpolationMethod::Linear => {
            let upper_weight = (value - values[lower]) / (values[upper] - values[lower]);
            stencil.push(lower, 1.0 - upper_weight);
            stencil.push(upper, upper_weight);
        }
        InterpolationMethod::Shell => {
            stencil.push(lower, 0.5);
            stencil.push(upper, 0.5);
        }
        InterpolationMethod::Lower => stencil.push(lower, 1.0),
    }
    stencil
}

fn unwrap_angle_near(mut angle: f64, center: f64) -> f64 {
    while angle - center > std::f64::consts::PI {
        angle -= 2.0 * std::f64::consts::PI;
    }
    while angle - center < -std::f64::consts::PI {
        angle += 2.0 * std::f64::consts::PI;
    }
    angle
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::raytracer::grid::{InterpolationWeight, MediumGrid};
    use crate::raytracer::primitive::PrimitiveKind;

    fn grid() -> StructuredGrid2D {
        StructuredGrid2D::new(
            10.0,
            vec![0.0, 10.0, 20.0],
            vec![-0.5, 0.0, 0.5],
            InterpolationMethod::Linear,
        )
        .unwrap()
    }

    fn point(radius: f64, angle: f64) -> Vec3 {
        radius * AngularBasis::CANONICAL.radial_direction(angle)
    }

    fn assert_weights(actual: &[InterpolationWeight], expected: &[(usize, f64)]) {
        assert_eq!(actual.len(), expected.len());
        for (actual, (expected_index, expected_weight)) in actual.iter().zip(expected) {
            assert_eq!(actual.index, *expected_index);
            assert!((actual.weight - expected_weight).abs() < 1e-12);
        }
    }

    #[test]
    fn rejects_invalid_axes_and_radii() {
        assert_eq!(
            StructuredGrid2D::new(
                -1.0,
                vec![0.0, 10.0],
                vec![-0.5, 0.5],
                InterpolationMethod::Linear,
            ),
            Err(GridError::InvalidEarthRadius)
        );
        assert_eq!(
            StructuredGrid2D::new(
                10.0,
                vec![0.0],
                vec![-0.5, 0.5],
                InterpolationMethod::Linear,
            ),
            Err(GridError::TooFewAltitudes)
        );
        assert_eq!(
            StructuredGrid2D::new(
                10.0,
                vec![0.0, 10.0],
                vec![0.0],
                InterpolationMethod::Linear,
            ),
            Err(GridError::TooFewHorizontalAngles)
        );
        assert_eq!(
            StructuredGrid2D::new(
                10.0,
                vec![0.0, 10.0],
                vec![0.0, -0.1],
                InterpolationMethod::Linear,
            ),
            Err(GridError::NonIncreasingHorizontalAngle { index: 1 })
        );
        assert_eq!(
            StructuredGrid2D::new(
                10.0,
                vec![0.0, 10.0],
                vec![0.0, std::f64::consts::PI],
                InterpolationMethod::Linear,
            ),
            Err(GridError::HorizontalSpanTooWide)
        );
        assert_eq!(
            StructuredGrid2D::new(
                10.0,
                vec![0.0, f64::NAN],
                vec![-0.5, 0.5],
                InterpolationMethod::Linear,
            ),
            Err(GridError::NonFiniteAltitude { index: 1 })
        );
        assert_eq!(
            StructuredGrid2D::new(
                10.0,
                vec![-10.0, 0.0],
                vec![-0.5, 0.5],
                InterpolationMethod::Linear,
            ),
            Err(GridError::InvalidRadius { index: 0 })
        );
        assert_eq!(
            StructuredGrid2D::new(
                f64::MAX,
                vec![0.0, f64::MAX],
                vec![-0.5, 0.5],
                InterpolationMethod::Linear,
            ),
            Err(GridError::InvalidRadius { index: 1 })
        );
        assert_eq!(
            StructuredGrid2D::new(
                10.0,
                vec![0.0, 10.0],
                vec![-0.5, f64::INFINITY],
                InterpolationMethod::Linear,
            ),
            Err(GridError::NonFiniteHorizontalAngle { index: 1 })
        );
    }

    #[test]
    fn angular_basis_rejects_degenerate_vectors() {
        assert_eq!(
            AngularBasis::new(Vec3::UNIT_Z, Vec3::UNIT_Z),
            Err(GridError::InvalidAngularBasis)
        );
        assert_eq!(
            AngularBasis::new(Vec3::ZERO, Vec3::UNIT_X),
            Err(GridError::InvalidAngularBasis)
        );
        assert_eq!(
            AngularBasis::new(Vec3::new(f64::NAN, 0.0, 1.0), Vec3::UNIT_X),
            Err(GridError::InvalidAngularBasis)
        );
        assert_eq!(
            AngularBasis::new(Vec3::new(f64::MAX, 0.0, 0.0), Vec3::UNIT_Z),
            Err(GridError::InvalidAngularBasis)
        );
    }

    #[test]
    fn uses_altitude_fastest_flattening() {
        let grid = grid();

        assert_eq!(grid.num_locations(), 9);
        assert_eq!(grid.location_shape(), (3, 3));
        assert_eq!(grid.num_cells(), 4);
        assert_eq!(grid.location_index(0, 0), 0);
        assert_eq!(grid.location_index(2, 0), 2);
        assert_eq!(grid.location_index(0, 1), 3);
        assert_eq!(grid.location_index(2, 2), 8);
    }

    #[test]
    fn locates_altitude_and_horizontal_cells() {
        let grid = grid();

        assert_eq!(grid.locate_indices(point(15.0, -0.25), 1e-12), Some((0, 0)));
        assert_eq!(grid.locate_indices(point(25.0, 0.25), 1e-12), Some((1, 1)));
        assert_eq!(
            grid.locate_cell(point(25.0, 0.25)),
            Some(CellId::Structured2D {
                altitude_index: 1,
                horizontal_index: 1,
            })
        );
    }

    #[test]
    fn rejects_points_outside_either_axis() {
        let grid = grid();

        assert!(grid.locate_cell(point(9.0, 0.0)).is_none());
        assert!(grid.locate_cell(point(31.0, 0.0)).is_none());
        assert!(grid.locate_cell(point(15.0, -0.6)).is_none());
        assert!(grid.locate_cell(point(15.0, 0.6)).is_none());
        assert!(grid.interpolation_weights_at(point(15.0, 0.6)).is_empty());
    }

    #[test]
    fn builds_bilinear_weights() {
        let weights: Vec<_> = grid()
            .interpolation_weights_at(point(15.0, -0.25))
            .iter()
            .collect();

        assert_weights(&weights, &[(0, 0.25), (1, 0.25), (3, 0.25), (4, 0.25)]);
    }

    #[test]
    fn collapses_weights_on_exact_grid_point() {
        let weights: Vec<_> = grid()
            .interpolation_weights_at(point(20.0, 0.0))
            .iter()
            .collect();

        assert_eq!(
            weights,
            vec![InterpolationWeight {
                index: 4,
                weight: 1.0
            }]
        );
    }

    #[test]
    fn clamps_boundary_roundoff() {
        let weights: Vec<_> = grid()
            .interpolation_weights_at(point(10.0 - 5e-9, -0.5 - 5e-9))
            .iter()
            .collect();

        assert_eq!(
            weights,
            vec![InterpolationWeight {
                index: 0,
                weight: 1.0,
            }]
        );
    }

    #[test]
    fn combines_shell_and_horizontal_weights() {
        let grid = StructuredGrid2D::new(
            10.0,
            vec![0.0, 10.0, 20.0],
            vec![-0.5, 0.0, 0.5],
            InterpolationMethod::Shell,
        )
        .unwrap();
        let weights: Vec<_> = grid
            .interpolation_weights_at(point(12.5, -0.25))
            .iter()
            .collect();

        assert_eq!(weights.len(), 4);
        assert!(
            weights
                .iter()
                .all(|weight| (weight.weight - 0.25).abs() < 1e-12)
        );
    }

    #[test]
    fn combines_lower_and_horizontal_weights() {
        let grid = StructuredGrid2D::new(
            10.0,
            vec![0.0, 10.0, 20.0],
            vec![-0.5, 0.0, 0.5],
            InterpolationMethod::Lower,
        )
        .unwrap();
        let weights: Vec<_> = grid
            .interpolation_weights_at(point(15.0, -0.25))
            .iter()
            .collect();

        assert_weights(&weights, &[(0, 0.5), (3, 0.5)]);
    }

    #[test]
    fn unwraps_angles_near_pi() {
        let grid = StructuredGrid2D::new(
            10.0,
            vec![0.0, 10.0],
            vec![3.0, 3.1, 3.2],
            InterpolationMethod::Linear,
        )
        .unwrap();
        let angle = -3.1;

        assert!(
            (grid.horizontal_angle_at(point(15.0, angle)) - (angle + 2.0 * std::f64::consts::PI))
                .abs()
                < 1e-12
        );
        assert!(grid.locate_cell(point(15.0, angle)).is_some());
    }

    #[test]
    fn respects_rotated_basis() {
        let basis = AngularBasis::new(Vec3::UNIT_X, -Vec3::UNIT_Z).unwrap();
        let grid = StructuredGrid2D::new_with_basis(
            10.0,
            vec![0.0, 10.0],
            vec![-0.5, 0.0, 0.5],
            InterpolationMethod::Linear,
            basis,
        )
        .unwrap();
        let point = 15.0 * basis.radial_direction(0.25);

        assert!((grid.horizontal_angle_at(point) - 0.25).abs() < 1e-12);
        assert_eq!(grid.locate_indices(point, 1e-12), Some((0, 1)));
    }

    #[test]
    fn builds_spherical_and_angular_primitives() {
        let primitives = grid().primitives();

        assert_eq!(primitives.len(), 6);
        assert_eq!(primitives[0].kind(), PrimitiveKind::Sphere);
        assert!(primitives[0].boundary().is_surface());
        assert!(primitives[2].boundary().is_top_of_atmosphere());
        assert_eq!(primitives[3].kind(), PrimitiveKind::AngularPlane);
        assert_eq!(
            primitives[3].boundary(),
            BoundaryTag::Horizontal { index: 0 }
        );
        assert_eq!(
            primitives[5].boundary(),
            BoundaryTag::Horizontal { index: 2 }
        );
    }
}
