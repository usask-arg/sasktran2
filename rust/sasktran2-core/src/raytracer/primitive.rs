use super::grid::BoundaryTag;
use super::ray::Ray;
use super::vec3::Vec3;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PrimitiveId(pub usize);

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Crossing {
    Entering,
    Exiting,
    Tangent,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PrimitiveKind {
    Sphere,
    Plane,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Intersection {
    pub distance: f64,
    pub point: Vec3,
    pub normal: Vec3,
    pub primitive_id: PrimitiveId,
    pub primitive_kind: PrimitiveKind,
    pub boundary: BoundaryTag,
    pub crossing: Crossing,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Primitive {
    Sphere {
        id: PrimitiveId,
        center: Vec3,
        radius: f64,
        radius_squared: f64,
        boundary: BoundaryTag,
    },
    Plane {
        id: PrimitiveId,
        point: Vec3,
        normal: Vec3,
        boundary: BoundaryTag,
    },
}

#[derive(Debug, Clone, Copy)]
struct PrimitiveMeta {
    id: PrimitiveId,
    kind: PrimitiveKind,
    boundary: BoundaryTag,
}

impl Primitive {
    #[inline(always)]
    pub fn sphere(id: PrimitiveId, center: Vec3, radius: f64, boundary: BoundaryTag) -> Self {
        Self::Sphere {
            id,
            center,
            radius,
            radius_squared: radius * radius,
            boundary,
        }
    }

    #[inline(always)]
    pub fn plane(id: PrimitiveId, point: Vec3, normal: Vec3, boundary: BoundaryTag) -> Self {
        Self::Plane {
            id,
            point,
            normal: normal.normalized(),
            boundary,
        }
    }

    #[inline(always)]
    pub fn id(&self) -> PrimitiveId {
        match self {
            Self::Sphere { id, .. } | Self::Plane { id, .. } => *id,
        }
    }

    #[inline(always)]
    pub fn boundary(&self) -> BoundaryTag {
        match self {
            Self::Sphere { boundary, .. } | Self::Plane { boundary, .. } => *boundary,
        }
    }

    #[inline(always)]
    pub fn kind(&self) -> PrimitiveKind {
        match self {
            Self::Sphere { .. } => PrimitiveKind::Sphere,
            Self::Plane { .. } => PrimitiveKind::Plane,
        }
    }

    pub fn intersections(&self, ray: Ray, epsilon: f64, output: &mut Vec<Intersection>) {
        match *self {
            Self::Sphere {
                id,
                center,
                radius,
                radius_squared,
                boundary,
            } => intersect_sphere(
                ray,
                center,
                radius,
                radius_squared,
                PrimitiveMeta {
                    id,
                    kind: PrimitiveKind::Sphere,
                    boundary,
                },
                epsilon,
                output,
            ),
            Self::Plane {
                id,
                point,
                normal,
                boundary,
            } => intersect_plane(
                ray,
                point,
                normal,
                PrimitiveMeta {
                    id,
                    kind: PrimitiveKind::Plane,
                    boundary,
                },
                epsilon,
                output,
            ),
        }
    }
}

fn push_intersection(
    distance: f64,
    ray: Ray,
    normal: Vec3,
    meta: PrimitiveMeta,
    crossing: Crossing,
    epsilon: f64,
    output: &mut Vec<Intersection>,
) {
    if distance < -epsilon {
        return;
    }
    let distance = if distance.abs() <= epsilon {
        0.0
    } else {
        distance
    };
    output.push(Intersection {
        distance,
        point: ray.point_at(distance),
        normal,
        primitive_id: meta.id,
        primitive_kind: meta.kind,
        boundary: meta.boundary,
        crossing,
    });
}

fn sphere_crossing(normal: Vec3, ray: Ray, discriminant_sqrt: f64, epsilon: f64) -> Crossing {
    if discriminant_sqrt <= epsilon {
        Crossing::Tangent
    } else if normal.dot(ray.direction) < 0.0 {
        Crossing::Entering
    } else {
        Crossing::Exiting
    }
}

fn intersect_sphere(
    ray: Ray,
    center: Vec3,
    radius: f64,
    radius_squared: f64,
    meta: PrimitiveMeta,
    epsilon: f64,
    output: &mut Vec<Intersection>,
) {
    let oc = ray.origin - center;
    let half_b = oc.dot(ray.direction);
    let c = oc.norm_squared() - radius_squared;
    let discriminant = half_b * half_b - c;

    if discriminant < -epsilon {
        return;
    }

    let discriminant_sqrt = discriminant.max(0.0).sqrt();
    let t_near = -half_b - discriminant_sqrt;
    let t_far = -half_b + discriminant_sqrt;

    let p_near = ray.point_at(t_near.max(0.0));
    let n_near = (p_near - center) / radius;
    push_intersection(
        t_near,
        ray,
        n_near,
        meta,
        sphere_crossing(n_near, ray, discriminant_sqrt, epsilon),
        epsilon,
        output,
    );

    if (t_far - t_near).abs() > epsilon {
        let p_far = ray.point_at(t_far.max(0.0));
        let n_far = (p_far - center) / radius;
        push_intersection(
            t_far,
            ray,
            n_far,
            meta,
            sphere_crossing(n_far, ray, discriminant_sqrt, epsilon),
            epsilon,
            output,
        );
    }
}

fn intersect_plane(
    ray: Ray,
    point: Vec3,
    normal: Vec3,
    meta: PrimitiveMeta,
    epsilon: f64,
    output: &mut Vec<Intersection>,
) {
    let denom = normal.dot(ray.direction);
    if denom.abs() <= epsilon {
        return;
    }

    let distance = (point - ray.origin).dot(normal) / denom;
    let crossing = if denom < 0.0 {
        Crossing::Entering
    } else {
        Crossing::Exiting
    };
    push_intersection(distance, ray, normal, meta, crossing, epsilon, output);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn surface() -> BoundaryTag {
        BoundaryTag::Surface { altitude_index: 0 }
    }

    #[test]
    fn sphere_intersection_returns_near_and_far_hits() {
        let primitive = Primitive::sphere(PrimitiveId(0), Vec3::ZERO, 1.0, surface());
        let ray = Ray::new(Vec3::new(0.0, 0.0, 3.0), Vec3::new(0.0, 0.0, -1.0));
        let mut hits = Vec::new();

        primitive.intersections(ray, 1e-12, &mut hits);

        assert_eq!(hits.len(), 2);
        assert!((hits[0].distance - 2.0).abs() < 1e-12);
        assert!((hits[1].distance - 4.0).abs() < 1e-12);
        assert_eq!(hits[0].crossing, Crossing::Entering);
        assert_eq!(hits[1].crossing, Crossing::Exiting);
    }

    #[test]
    fn tangent_sphere_intersection_is_single_hit() {
        let primitive = Primitive::sphere(PrimitiveId(0), Vec3::ZERO, 1.0, surface());
        let ray = Ray::new(Vec3::new(1.0, 0.0, 3.0), Vec3::new(0.0, 0.0, -1.0));
        let mut hits = Vec::new();

        primitive.intersections(ray, 1e-12, &mut hits);

        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].crossing, Crossing::Tangent);
        assert!((hits[0].distance - 3.0).abs() < 1e-12);
    }

    #[test]
    fn plane_intersection_returns_forward_hit() {
        let primitive = Primitive::plane(PrimitiveId(0), Vec3::ZERO, Vec3::UNIT_Z, surface());
        let ray = Ray::new(Vec3::new(0.0, 0.0, 3.0), Vec3::new(0.0, 0.0, -1.0));
        let mut hits = Vec::new();

        primitive.intersections(ray, 1e-12, &mut hits);

        assert_eq!(hits.len(), 1);
        assert!((hits[0].distance - 3.0).abs() < 1e-12);
    }
}
