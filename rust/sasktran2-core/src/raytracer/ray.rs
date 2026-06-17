use super::primitive::{Intersection, Primitive};
use super::vec3::Vec3;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ray {
    pub origin: Vec3,
    pub direction: Vec3,
}

impl Ray {
    #[inline(always)]
    pub fn new(origin: Vec3, direction: Vec3) -> Self {
        Self {
            origin,
            direction: direction.normalized(),
        }
    }

    #[inline(always)]
    pub fn point_at(self, distance: f64) -> Vec3 {
        self.origin + self.direction * distance
    }

    #[inline(always)]
    pub fn cos_zenith_at_origin(self) -> f64 {
        self.origin.normalized().dot(self.direction)
    }

    #[inline(always)]
    pub fn spherical_tangent_distance(self) -> f64 {
        -self.origin.dot(self.direction)
    }

    #[inline(always)]
    pub fn spherical_tangent_radius(self) -> f64 {
        self.origin.cross(self.direction).norm()
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StraightPath {
    ray: Ray,
}

impl StraightPath {
    #[inline(always)]
    pub fn new(ray: Ray) -> Self {
        Self { ray }
    }

    #[inline(always)]
    pub fn position(&self, distance: f64) -> Vec3 {
        self.ray.point_at(distance)
    }

    #[inline(always)]
    pub fn direction(&self, _distance: f64) -> Vec3 {
        self.ray.direction
    }

    #[inline(always)]
    pub fn intersections(
        &self,
        primitive: &Primitive,
        epsilon: f64,
        output: &mut Vec<Intersection>,
    ) {
        primitive.intersections(self.ray, epsilon, output);
    }
}
